#define _POSIX_C_SOURCE 200809L
#include <sys/types.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* C (n x p) = A (n x m) * B (m x p), on a q x q periodic process grid.
 * Every dimension is split into q nearly-equal pieces (first dim%q pieces get +1).
 * Each local block lives in a zero-padded, fixed-size buffer so that all
 * messages have a constant size. */

static inline int blk_size(int len, int q, int t) { return len / q + (t < len % q); }
static inline int blk_off (int len, int q, int t) { return t * (len / q) + (t < len % q ? t : len % q); }
static inline int blk_max (int len, int q)        { return (len + q - 1) / q; }

/* ---------- Matrix reading from file ---------- */
void read_matrix_streaming(const char *filename, double **A, double **B, int *n, int *m, int *p) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("fopen"); MPI_Abort(MPI_COMM_WORLD, 1); }

    char *line = NULL;
    size_t len = 0;
    int readingA = 0, readingB = 0;
    int rowsA = 0, colsA = 0, rowsB = 0, colsB = 0;

    /* First pass: dimensions */
    while (getline(&line, &len, fp) != -1) {
        if (strncmp(line, "Matrix A", 8) == 0) { readingA = 1; readingB = 0; continue; }
        if (strncmp(line, "Matrix B", 8) == 0) { readingA = 0; readingB = 1; continue; }
        if (!readingA && !readingB) continue;

        int count = 0;
        char *tmp = strtok(line, " \t\n");
        while (tmp) { count++; tmp = strtok(NULL, " \t\n"); }
        if (count == 0) continue;

        if (readingA) { rowsA++; colsA = count; }
        else          { rowsB++; colsB = count; }
    }

    if (colsA != rowsB) {
        fprintf(stderr, "Dimension mismatch: A is %dx%d, B is %dx%d\n", rowsA, colsA, rowsB, colsB);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    rewind(fp);

    *n = rowsA;
    *m = colsA;
    *p = colsB;
    *A = malloc((size_t)rowsA * colsA * sizeof(double));
    *B = malloc((size_t)rowsB * colsB * sizeof(double));

    size_t idxA = 0, idxB = 0;
    readingA = readingB = 0;

    /* Second pass: values */
    while (getline(&line, &len, fp) != -1) {
        if (strncmp(line, "Matrix A", 8) == 0) { readingA = 1; readingB = 0; continue; }
        if (strncmp(line, "Matrix B", 8) == 0) { readingA = 0; readingB = 1; continue; }
        if (!readingA && !readingB) continue;

        char *tmp = strtok(line, " \t\n");
        while (tmp) {
            double val = atof(tmp);
            if (readingA) (*A)[idxA++] = val;
            else          (*B)[idxB++] = val;
            tmp = strtok(NULL, " \t\n");
        }
    }

    free(line);
    fclose(fp);
}

/* ---------- Local multiply on padded buffers: C(mi x nj) += A(mi x kc) * B(kc x nj) ---------- */
static void local_multiply(int mi, int nj, int kc,
                           const double *A, int lda,
                           const double *B, int ldb,
                           double *C, int ldc)
{
    for (int i = 0; i < mi; i++)
        for (int k = 0; k < kc; k++) {
            double aik = A[(size_t)i * lda + k];
            for (int j = 0; j < nj; j++)
                C[(size_t)i * ldc + j] += aik * B[(size_t)k * ldb + j];
        }
}

/* ---------- Write matrix ---------- */
void write_matrix_to_file(const char *filename, const double *C, int n, int p) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            fprintf(fp, "%.4f", C[(size_t)i * p + j]);
            if (j < p - 1) fprintf(fp, " ");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void PrintMatrix(int n, int m, const double *mat) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            printf("%8.2f ", mat[(size_t)i * m + j]);
        printf("\n");
    }
}

/* ---------- Square Cartesian grid ---------- */
MPI_Comm create_cartesian_grid(MPI_Comm base_comm, int *q_out)
{
    int npes, rank;
    MPI_Comm_size(base_comm, &npes);
    MPI_Comm_rank(base_comm, &rank);

    int q = (int)(sqrt((double)npes) + 0.5);
    if (q * q != npes) {
        if (rank == 0)
            fprintf(stderr, "Number of processes (%d) must be a perfect square\n", npes);
        MPI_Abort(base_comm, 1);
    }

    int dims[2]    = {q, q};
    int periods[2] = {1, 1};
    MPI_Comm grid;
    /* reorder = 0 so grid rank 0 == world rank 0 (the one that read the file) */
    MPI_Cart_create(base_comm, 2, dims, periods, 0, &grid);

    *q_out = q;
    return grid;
}

/* ---------- Scatter a rows x cols matrix into padded (rbm x cbm) blocks ---------- */
void scatter_blocks(const double *G, int rows, int cols,
                    double *local, int rbm, int cbm, int q, MPI_Comm grid)
{
    int rank, npes;
    MPI_Comm_rank(grid, &rank);
    MPI_Comm_size(grid, &npes);
    size_t bufsz = (size_t)rbm * cbm;

    if (rank == 0) {
        double *buf = malloc(bufsz * sizeof(double));
        for (int r = 0; r < npes; r++) {
            int rc[2];
            MPI_Cart_coords(grid, r, 2, rc);
            int br = blk_size(rows, q, rc[0]), bc = blk_size(cols, q, rc[1]);
            int ro = blk_off (rows, q, rc[0]), co = blk_off (cols, q, rc[1]);

            double *dst = (r == 0) ? local : buf;
            memset(dst, 0, bufsz * sizeof(double));
            for (int i = 0; i < br; i++)
                memcpy(dst + (size_t)i * cbm, G + (size_t)(ro + i) * cols + co, bc * sizeof(double));

            if (r != 0)
                MPI_Send(buf, (int)bufsz, MPI_DOUBLE, r, 999, grid);
        }
        free(buf);
    } else {
        MPI_Recv(local, (int)bufsz, MPI_DOUBLE, 0, 999, grid, MPI_STATUS_IGNORE);
    }
}

/* ---------- Gather padded blocks into a rows x cols matrix on rank 0 ---------- */
void gather_blocks(const double *local, double *G, int rows, int cols,
                   int rbm, int cbm, int q, MPI_Comm grid)
{
    int rank, npes;
    MPI_Comm_rank(grid, &rank);
    MPI_Comm_size(grid, &npes);
    size_t bufsz = (size_t)rbm * cbm;

    if (rank == 0) {
        double *buf = malloc(bufsz * sizeof(double));
        for (int r = 0; r < npes; r++) {
            int rc[2];
            MPI_Cart_coords(grid, r, 2, rc);
            int br = blk_size(rows, q, rc[0]), bc = blk_size(cols, q, rc[1]);
            int ro = blk_off (rows, q, rc[0]), co = blk_off (cols, q, rc[1]);

            const double *src = local;
            if (r != 0) {
                MPI_Recv(buf, (int)bufsz, MPI_DOUBLE, r, 1999, grid, MPI_STATUS_IGNORE);
                src = buf;
            }
            for (int i = 0; i < br; i++)
                memcpy(G + (size_t)(ro + i) * cols + co, src + (size_t)i * cbm, bc * sizeof(double));
        }
        free(buf);
    } else {
        MPI_Send(local, (int)bufsz, MPI_DOUBLE, 0, 1999, grid);
    }
}

/* ---------- Cannon: C (n x p) = A (n x m) * B (m x p) ---------- */
void MatrixMatrixMultiply(int n, int m, int p,
                          double *a, double *b, double *c, MPI_Comm grid)
{
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(grid, 2, dims, periods, coords);
    int q = dims[0];
    int i = coords[0], j = coords[1];

    int nb = blk_max(n, q), mb = blk_max(m, q), pb = blk_max(p, q);  /* padded block dims */
    int ni = blk_size(n, q, i);   /* my rows of A and C   */
    int pj = blk_size(p, q, j);   /* my cols of B and C   */

    int left, right, up, down, src, dst;
    MPI_Cart_shift(grid, 1, +1, &left, &right);   /* A moves left  */
    MPI_Cart_shift(grid, 0, +1, &up,   &down);    /* B moves up    */

    /* Initial skew: A row i left by i, B column j up by j */
    MPI_Cart_shift(grid, 1, -i, &src, &dst);
    MPI_Sendrecv_replace(a, nb * mb, MPI_DOUBLE, dst, 100, src, 100, grid, MPI_STATUS_IGNORE);
    MPI_Cart_shift(grid, 0, -j, &src, &dst);
    MPI_Sendrecv_replace(b, mb * pb, MPI_DOUBLE, dst, 200, src, 200, grid, MPI_STATUS_IGNORE);

    for (int s = 0; s < q; s++) {
        int kidx = (i + j + s) % q;          /* k-block currently held (both A and B) */
        int kc   = blk_size(m, q, kidx);     /* its true extent */

        local_multiply(ni, pj, kc, a, mb, b, pb, c, pb);

        MPI_Sendrecv_replace(a, nb * mb, MPI_DOUBLE, left, 1000 + s, right, 1000 + s, grid, MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(b, mb * pb, MPI_DOUBLE, up,   2000 + s, down,  2000 + s, grid, MPI_STATUS_IGNORE);
    }
    /* No need to restore the original distribution unless you reuse a and b. */
}

/* ---------- Main ---------- */
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, npes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    if (argc < 3) {
        if (rank == 0)
            fprintf(stderr, "Usage: %s <input_file> <output_file>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    const char *input_file  = argv[1];
    const char *output_file = argv[2];

    double *A = NULL, *B = NULL, *C_full = NULL;
    int dimsnmp[3] = {0, 0, 0};   /* n, m, p */

    if (rank == 0)
        read_matrix_streaming(input_file, &A, &B, &dimsnmp[0], &dimsnmp[1], &dimsnmp[2]);

    MPI_Bcast(dimsnmp, 3, MPI_INT, 0, MPI_COMM_WORLD);
    int n = dimsnmp[0], m = dimsnmp[1], p = dimsnmp[2];

    int q;
    MPI_Comm grid = create_cartesian_grid(MPI_COMM_WORLD, &q);

    int nb = blk_max(n, q), mb = blk_max(m, q), pb = blk_max(p, q);

    double *a_local = calloc((size_t)nb * mb, sizeof(double));
    double *b_local = calloc((size_t)mb * pb, sizeof(double));
    double *c_local = calloc((size_t)nb * pb, sizeof(double));

    scatter_blocks(A, n, m, a_local, nb, mb, q, grid);
    scatter_blocks(B, m, p, b_local, mb, pb, q, grid);

    MPI_Barrier(grid);
    double start = MPI_Wtime();
    MatrixMatrixMultiply(n, m, p, a_local, b_local, c_local, grid);
    MPI_Barrier(grid);
    double end = MPI_Wtime();
    if (rank == 0)
        printf("Time for matrix multiplication: %.6f seconds\n", end - start);

    if (rank == 0)
        C_full = malloc((size_t)n * p * sizeof(double));
    gather_blocks(c_local, C_full, n, p, nb, pb, q, grid);

    if (rank == 0) {
        if (n <= 16 && m <= 16 && p <= 16) {   /* only dump small matrices */
            printf("A\n"); PrintMatrix(n, m, A);
            printf("B\n"); PrintMatrix(m, p, B);
            printf("C\n"); PrintMatrix(n, p, C_full);
        }
        write_matrix_to_file(output_file, C_full, n, p);
        printf("Result written to %s\n", output_file);
        free(A);
        free(B);
        free(C_full);
    }

    MPI_Comm_free(&grid);
    free(a_local);
    free(b_local);
    free(c_local);

    MPI_Finalize();
    return 0;
}
