#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


void get_code_id(char* code_id)
{
   //string is replaced by script generate_userblock.sh
   strcpy(code_id,"DUMMY_INSERT_CODE_ID_HERE");
}

void get_build_id(char* build_id)
{
   //string is replaced by script generate_userblock.sh
   strcpy(build_id,"DUMMY_INSERT_BUILD_ID_HERE");
}


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


void check_file_exists(char* subdirname, char* filename, bool* exists)
{
   struct stat st = {0};
   if (stat(subdirname, &st) == -1) {
      mkdir(subdirname, 0775);
   }
   char filepath[255];
   strcpy(filepath,subdirname);
   strcat(filepath,"/");
   strcat(filepath,filename);
   *exists = stat(filepath, &st) != -1;
}


void print_userblock(char* filename, char* inifilename)
{
   FILE* fp = fopen(filename, "w");
   rewind(fp);

   fprintf(fp, "{[( START USERBLOCK )]}\n");

   fprintf(fp, "{[( INIFILE )]}\n");
   struct stat st = {0};
   if ( stat(inifilename, &st) != -1 ) {
      FILE* fini = fopen(inifilename, "rb");
      int c;
      do {
         c = fgetc (fini);
         if (c != EOF) fputc((char)c, fp);
      } while (c != EOF);
      fclose(fini);
   }

   /* the shell script generate_userblock.sh copies
      this file into the build directory and inserts
      print commands for the code and build info here. */

   //INSERT_CODE_INFO_HERE

   //INSERT_BUILD_INFO_HERE

   fprintf(fp, "{[( END USERBLOCK )]}\n");

   fclose(fp);
}


void print_run_info(char* filename, char* inifilename, char* runid, char* command)
{
   FILE* fp = fopen(filename, "w");
   rewind(fp);

   fprintf(fp, "{[( RUN ID )]}\n");
   fprintf(fp, "%s\n", runid);

   fprintf(fp, "{[( BUILD ID )]}\n");
   //string is replaced by script generate_userblock.sh
   fprintf(fp, "DUMMY_INSERT_BUILD_ID_HERE\n");

   fprintf(fp, "{[( CODE ID )]}\n");
   //string is replaced by script generate_userblock.sh
   fprintf(fp, "DUMMY_INSERT_CODE_ID_HERE\n");

   fprintf(fp, "{[( COMMAND )]}\n");
   fprintf(fp, "%s\n", command);
   fprintf(fp, "{[( INIFILE )]}\n");

   struct stat st = {0};
   if ( stat(inifilename, &st) != -1 ) {
      FILE* fini = fopen(inifilename, "rb");
      int c;
      do {
         c = fgetc (fini);
         if (c != EOF) fputc((char)c, fp);
      } while (c != EOF);
      fclose(fini);
   }

   fprintf(fp, "{[( DEPENDENCIES )]}\n");

   fclose(fp);
}


void print_code_info(char* filename)
{

   FILE* fp = fopen(filename, "w");
   rewind(fp);

   fprintf(fp, "{[( CODE ID )]}\n");
   //string is replaced by script generate_userblock.sh
   fprintf(fp, "DUMMY_INSERT_CODE_ID_HERE\n");

   /* the shell script generate_userblock.sh copies
      this file into the build directory and inserts
      print commands for the code info here. */

   //INSERT_CODE_INFO_HERE

   fclose(fp);
}


void print_build_info(char* filename)
{

   FILE* fp = fopen(filename, "w");
   rewind(fp);

   fprintf(fp, "{[( BUILD ID )]}\n");
   //string is replaced by script generate_userblock.sh
   fprintf(fp, "DUMMY_INSERT_BUILD_ID_HERE\n");

   fprintf(fp, "{[( CODE ID )]}\n");
   //string is replaced by script generate_userblock.sh
   fprintf(fp, "DUMMY_INSERT_CODE_ID_HERE\n");

   /* the shell script generate_userblock.sh copies
      this file into the build directory and inserts
      print commands for the build info here. */

   //INSERT_BUILD_INFO_HERE

   fclose(fp);
}
