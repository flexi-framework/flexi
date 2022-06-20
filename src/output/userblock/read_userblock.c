#include <string.h>
#include <stdio.h>

void copy_userblock(char* outfilename, char* infilename)
{
   FILE* fout = fopen(outfilename,"rb+");
   FILE* fin  = fopen(infilename, "r");
   rewind(fout);
   int c;
   do {
      c = fgetc (fin);
      if (c != EOF) fputc((char) c, fout);
   } while (c != EOF);
   fclose(fin);
   fclose(fout);
}

void print_userblock(char* filename, char* inifilename)
{
   FILE* fp = fopen(filename, "w");
   rewind(fp);
   fprintf(fp, "{[( START USERBLOCK )]}\n");

   // ini file
   fprintf(fp, "{[( INIFILE )]}\n");
   FILE* fini = fopen(inifilename, "rb");
   int c;
   do {
      c = fgetc (fini);
      if (c != EOF) fputc((char)c, fp);
   } while (c != EOF);
   fclose(fini);

   /*the shell script generate_userblock.sh copies
    * this file into the build directory and inserts
    * print commands for the build info here. */

   //INSERT_BUILD_INFO_HERE

   fprintf(fp, "{[( END USERBLOCK )]}\n");
   fclose(fp);
}
