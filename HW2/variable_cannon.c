#define _POSIX_C_SOURCE 200809L  // enables getline() and ssize_t on most compilers
#include <sys/types.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* ---------- Matrix reading from file (use streaming for very large files >500MB) ---------- */
static inline int blk_size(int len, int q, int t) { return len / q + (t < len % q); }
static inline int blk_off (int len, int q, int t) { return t * (len / q) + (t < len % q ? t : len % q); }
static inline int blk_max (int len, int q)        { return (len + q - 1) / q; }
void read_matrix_streaming(const char *filename, double **A, double **B, int *n, int *m, int *p) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("fopen"); MPI_Abort(MPI_COMM_WORLD, 1); }

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int readingA = 0, readingB = 0;
    int rowsA = 0, colsA = 0, rowsB = 0, colsB = 0;

    // First pass — determine dimensions
    while ((read = getline(&line, &len, fp)) != -1) {
        if (strncmp(line, "Matrix A", 8) == 0) { readingA = 1; readingB = 0; continue; }
        if (strncmp(line, "Matrix B", 8) == 0) { readingA = 0; readingB = 1; continue; }
        if (!readingA && !readingB) continue;

        int count = 0;
        char *tmp = strtok(line, " \t\n");
        while (tmp) { count++; tmp = strtok(NULL, " \t\n"); }
        if (count == 0) continue;

        if (readingA) { rowsA++; colsA = count; }
        else if (readingB) { rowsB++; colsB = count; }
    }

    if(colsA != rowsB){
       fprintf(stderr, "cant do those dimensions bucko");
       MPI_Abort(MPI_COMM_WORLD, 1);
    }

    rewind(fp);

    *n = rowsA;
    *m = colsA;
    *p = colsB;
    *A = malloc((size_t)rowsA * colsA * sizeof(double));
    *B = malloc((size_t)rowsB * colsB * sizeof(double));

    int idxA = 0, idxB = 0;
    readingA = readingB = 0;

    // Second pass — parse data directly into flat arrays
    while ((read = getline(&line, &len, fp)) != -1) {
        if (strncmp(line, "Matrix A", 8) == 0) { readingA = 1; readingB = 0; continue; }
        if (strncmp(line, "Matrix B", 8) == 0) { readingA = 0; readingB = 1; continue; }
        if (!readingA && !readingB) continue;

        char *tmp = strtok(line, " \t\n");
        while (tmp) {
            double val = atof(tmp);
            if (readingA) (*A)[idxA++] = val;
            else if (readingB) (*B)[idxB++] = val;
            tmp = strtok(NULL, " \t\n");
        }
    }

    free(line);
    fclose(fp);
}


/* ---------- Local matrix multiplication ---------- */
void MatrixMultiply(int nl, int ml, int pl, double *A, int lda, double *B, int ldb,  double *C, int ldc) {
    for (int lv1 = 0; lv1 < nl; lv1++) {
        const double *Ai = A + (size_t)lv1 * lda;
        double *Ci = C + (size_t)lv1 * ldc;
        for (int lv3 = 0; lv3 < pl; lv3++) {
            double a = Ai[lv3];
            const double *Bk = B + (size_t)lv3 * ldb;
            for (int lv2 = 0; lv2 < ml; lv2++)
                Ci[lv2] += a * Bk[lv2];
        }
    }
}

/* ---------- Write matrix to output file ---------- */
void write_matrix_to_file(const char *filename, double *C, int n, int p) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            fprintf(fp, "%.4f", C[i * p + j]);
            if (j < p - 1) fprintf(fp, " ");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

/* ---------- Debug: Print matrix ---------- */
void PrintMatrix(int n, int m, double* mat) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%8.2f ", mat[i * m + j]);
        }
        printf("\n");
    }
}

/* ---------- Create Cartesian grid communicator ---------- */
MPI_Comm create_cartesian_grid(MPI_Comm base_comm, int *q_out)
{
    int npes, rank;
    MPI_Comm_size(base_comm, &npes);
    MPI_Comm_rank(base_comm, &rank);

    int dims[2] = {0, 0};
    MPI_Dims_create(npes, 2, dims);
    int q = (int)(sqrt((double)npes) + 0.5);
    if (q * q != npes) {
        if (rank == 0)
            fprintf(stderr, "Number of processes (%d) must be a perfect square\n", npes);
        MPI_Abort(base_comm, 1);
    }

    int periods[2] = {1, 1};  /* wraparound grid for Cannon */
    MPI_Comm grid;
    MPI_Cart_create(base_comm, 2, dims, periods, 0, &grid);

    *q_out = q;
    return grid;
}

/* ---------- Scatter global n×n matrix into nlocal×nlocal blocks ---------- */
void scatter_2d_blocks(double *A_global, double *local, int rows, int cols, int rbm, int cbm, int q, MPI_Comm grid){
    int npes, rank;
    MPI_Comm_size(grid, &npes);
    MPI_Comm_rank(grid, &rank);

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
                memcpy(dst + (size_t)i * cbm, A_global + (size_t)(ro + i) * cols + co, bc * sizeof(double));

            if (r != 0)
                MPI_Send(buf, (int)bufsz, MPI_DOUBLE, r, 999, grid);
        }
        free(buf);
    } else {
        MPI_Recv(local, (int)bufsz, MPI_DOUBLE, 0, 999, grid, MPI_STATUS_IGNORE);
    }
}

/* ---------- Gather nlocal×nlocal blocks back into global n×n matrix ---------- */
void gather_2d_blocks(double *local, double *C_global, int rows, int cols, int rbm, int cbm, int q, MPI_Comm grid){
    int npes, rank;
    MPI_Comm_size(grid, &npes);
    MPI_Comm_rank(grid, &rank);
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
                memcpy(C_global + (size_t)(ro + i) * cols + co, src + (size_t)i * cbm, bc * sizeof(double));
        }
        free(buf);
    } else {
        MPI_Send(local, (int)bufsz, MPI_DOUBLE, 0, 1999, grid);
    }
}

void MatrixMatrixMultiply(int n, int m, int p, double *a, double *b, double *c, MPI_Comm grid)
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

        MatrixMultiply(ni, pj, kc, a, mb, b, pb, c, pb);

        MPI_Sendrecv_replace(a, nb * mb, MPI_DOUBLE, left, 1000 + s, right, 1000 + s, grid, MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(b, mb * pb, MPI_DOUBLE, up,   2000 + s, down,  2000 + s, grid, MPI_STATUS_IGNORE);
    }
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

    const char *input_file = argv[1];
    const char *output_file = argv[2];

    double *A = NULL, *B = NULL, *C_full = NULL;
    int dimsnmp[3] = {0, 0, 0};

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

    scatter_2d_blocks(A, a_local, n, m, nb, mb, q, grid);
    scatter_2d_blocks(B, b_local, m, p, mb, pb, q, grid);

    double start, end, elapsed;
    start = MPI_Wtime();
    /* Run Cannon’s algorithm */
    MatrixMatrixMultiply(n, m, p, a_local, b_local, c_local, grid);
    end = MPI_Wtime();
    if (rank == 0) {
        elapsed = end - start;
        printf("Time for matrix multiplication: %.6f seconds\n", elapsed);
    }

    /* Gather result back */
    if (rank == 0)
        C_full = malloc((size_t)n * p * sizeof(double));
    gather_2d_blocks(c_local, C_full, n, p, nb, pb, q, grid);

    /*if (rank == 0){
        printf("A\n");
        PrintMatrix(n, m, A);
        printf("B\n");
        PrintMatrix(m, p, B);
        printf("C\n");
        PrintMatrix(n, p, C_full);
    }*/
    /* Write result to output file */
    if (rank == 0) {
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
