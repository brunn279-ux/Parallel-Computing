#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

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


int main(int argc, char *argv[]){


   if (argc < 2){

      fprintf(stderr, "include a filename");
      return 1;
   }

   FILE *inpf = fopen(argv[1], "r");
   if (inpf == NULL){

      perror("Error opening file");
      return 1;

   }


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

   printf("A is %dx%d, B is %dx%d\n", A.rows, A.cols, B.rows, B.cols);

   free(A.data);
   free(B.data);

   return 0;
}
