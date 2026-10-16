#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define LINE_MAX 8192

typedef struct{

   double *data;
   int rows, cols;

} Matrix;

static int is_blank(const char *s){

   while (*s){
      if(!isspace((unsigned char)*s++)){
         return 0;
      }
   }
   return 1;

}

static int is_label(const char *s){
   while(isspace((unsigned char)*s)){
      s++;
   }
   return !(isdigit((unsigned char)*s) || *s == '-' || *s == '+' || *s == '.');
}

static int count_values(const char *s){
   int n = 0;
   char *end;
   for (;;){
      strtod(s, &end);
      if (end == s) break;
      n++;
      s = end;
   }
   while(isspace((unsigned char)*s)) s++;
   return *s == '\0' ? n : -1;
}

int load_mats(FILE *inpf, Matrix *A, Matrix *B){

   Matrix *mats[2] = {A, B};
   size_t cap[2] = {0, 0};
   int count = 0;
   Matrix *cur = NULL;
   size_t *curcap = NULL;
   char line[LINE_MAX];

   A->data = B->data = NULL;
   A->rows = A->cols = B->rows = B->cols  = 0;

   while(fgets(line, sizeof line, inpf)){
      if(is_blank(line)) continue;

      if (is_label(line)){
         if(count >= 2){
            fprintf(stderr, "more than two mats\n");
            return -1;
         }
         cur = mats[count];
         curcap = &cap[count];
         count++;
         continue;
      }

      if (cur == NULL){
         fprintf(stderr, "data before label\n");
         return -1;
      }

      int n = count_values(line);
      if (n <= 0){

         fprintf(stderr, "bad data line: %s", line);
         return -1;


      }
      if (cur->cols == 0)
         cur->cols = n;
      else if (n != cur->cols){
         fprintf(stderr, "row %d has %d vals, expected %d\n", cur->rows + 1, n, cur->cols);
         return -1;
      }

      size_t need = (size_t)(cur->rows + 1) * cur->cols;
      if (need > *curcap){
         size_t newcap = *curcap ? *curcap * 2 : (size_t)cur->cols * 8;
         while (newcap < need) newcap *= 2;
         double *tmp = realloc(cur->data, newcap * sizeof *tmp);
         if (tmp==NULL){
            perror("realloc");
            return -1;
         }
         cur->data = tmp;
         *curcap = newcap;
      }

      char *p = line, *end;
      for (int lv1 = 0; lv1 < cur->cols; lv1++){
         cur->data[(size_t)cur->rows * cur->cols + lv1] = strtod(p, &end);
         p = end;
      }
      cur->rows++;

   }

   if(count != 2 || A->rows == 0 || B->rows == 0){
      fprintf(stderr, "expected two non-empty matricies, found %d\n", count);
      return -1;
   }

   return 0;


}


double norm(Matrix *A, Matrix *B, int aind, int bind){

   double sum = 0;
   double tempA, tempB;

   for(int lv1 = 0; lv1 < A->cols; lv1++){
      tempA = A -> data[aind * A->cols + lv1];
      tempB = B -> data[bind * B->cols + lv1];

      sum += (tempA - tempB)*(tempA - tempB);
   }

   return sqrt(sum);

}


double minVal(double a, double b){
   if(a<b){
      return a;
   }
   else{
      return b;
   }
}

double DTW(Matrix *A, Matrix *B, double *elapsed){

   int m = A->rows, p = B->rows;
   double (*D)[p] = calloc(m, sizeof *D);
   double (*dist)[p] = calloc(m, sizeof *dist);

   double cost = 0;

   if (D == NULL || dist == NULL){
      perror("calloc");
      free(D);
      free(dist);
      return -1;
   }

   double start_time = omp_get_wtime();

   #pragma omp parallel for collapse(2)
   for (int lv1 = 0; lv1 < A->rows; lv1++){
      for (int lv2 = 0; lv2 < B->rows; lv2++){
         dist[lv1][lv2] = norm(A, B, lv1, lv2);
      }
   }

   D[0][0] = dist[0][0];

   for (int lv1 = 1; lv1 < A->rows; lv1++){
      D[lv1][0] = dist[lv1][0]+ D[lv1-1][0];
   }
   for (int lv1 = 1; lv1 < B->rows; lv1++){
      D[0][lv1] = dist[0][lv1] + D[0][lv1-1];
   }
   /*for(int k = 1; k < 2*(arrSize)- 1; k++){
      if(k <= arrSize){
      #pragma omp parallel for collapse(2) private(lv1, lv2)
         for(int lv1 = 1; lv1 <= k; lv1++){
            for(int lv2 = k; lv2 >= 1; lv2--){
               D[lv1][lv2] = dist[lv1][lv2] + minVal(minVal(D[lv1-1][lv2], D[lv1][lv2-1]), D[lv1-1][lv2-1]);
            }
         }
      }
      else if(k > arrSize){
      #pragma omp parallel for collapse(2) private(lv1, lv2)
         for(int lv1 = arrSize - (2*arrSize-k); lv1 <= arrSize; lv1++){
            for(int lv2 = arrSize; lv2 >= arrSize - (2*arrSize-k); lv2--){
               D[lv1][lv2] = dist[lv1][lv2] + minVal(minVal(D[lv1-1][lv2], D[lv1][lv2-1]), D[lv1-1][lv2-1]);
            }
         }
      }
      else if (k > 2*arrSize -1){
         perror("exceded array bounds");
         return -1;
      }
   }
   */

   #pragma omp parallel

   for (int k = 2; k <= m + p - 2; k++){
      int lv1Min = (k-(p-1) > 1) ? k-(p-1) : 1;
      int lv1Max = (k-1 < m -1) ? k-1 : m-1;

      #pragma omp for
      for (int lv1 = lv1Min; lv1 <= lv1Max; lv1++){
         int lv2 = k-lv1;

         D[lv1][lv2] = dist[lv1][lv2] + minVal(minVal(D[lv1-1][lv2], D[lv1][lv2-1]), D[lv1-1][lv2-1]);
      }
   }

   double end_time = omp_get_wtime();

   cost = D[m-1][p-1];
   *elapsed = end_time - start_time;

   free(D);
   free(dist);

   return cost;
}

int main(int argc, char *argv[]){

   omp_set_dynamic(0);

   double alCost = 0;
   double tTime = 0;

   int tcount = 0;

   if (argc < 2){

      fprintf(stderr, "include a filename");
      return 1;
   }

   FILE *inpf = fopen(argv[1], "r");
   if (inpf == NULL){

      perror("Error opening file");
      return 1;

   }

   tcount = atoi(argv[2]);
   omp_set_num_threads(tcount);
   //printf("thread count is %d\n", tcount);

   Matrix A, B;
   int rc = load_mats(inpf, &A, &B);
   fclose(inpf);

   if (rc != 0){
      free(A.data);
      free(B.data);
      return 1;
   }
   if (A.cols != B.cols){
      fprintf(stderr, "column mismatch A:%d, B:%d\n", A.cols, B.cols);
      free(A.data);
      free(B.data);
      return 1;
   }

   //printf("A is %dx%d, B is %dx%d\n", A.rows, A.cols, B.rows, B.cols);

   alCost = DTW(&A, &B, &tTime);

   printf("The Final Cost is: %f\n", alCost);
   printf("The Total Completion time (ms) is: %f\n", tTime*1000);

   free(A.data);
   free(B.data);

   return 0;
}
