/*=================================================================================================================================
 * Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
 * This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
 * For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
 *
 * FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
 *
 * You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
 *================================================================================================================================*/


/* ====================================================================
 * The shell script "generate_userblock.sh" copies this file into the
 * build directory and inserts print commands into this source file
 * specific information about the code version and build configuration.
 * This contains specifically:
 *   - Code information (Git branch, commit and diff)
 *   - Build information (CMake configuration)
 *
 * Moreover, this file implements functionality regarding the
 * "userblock" of FLEXI, which appends this information to each result
 * file of FLEXI, such that each file contains the code and build
 * configuration as well as the used parameter file, in order to
 * reproduce the obtained results.
 * ===================================================================*/


#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


/* Returns the unique code ID as string, which is determined and
 * inserted here during compile time by the "generate_userblock.sh"
 * script.
 */
void get_code_id(char* code_id)
{
   // Dummy string is replaced by script "generate_userblock.sh"
   strcpy(code_id,"DUMMY_INSERT_CODE_ID_HERE");
}


/* Returns the unique build ID as string, which is determined and
 * inserted here during compile time by the "generate_userblock.sh"
 * script.
 */
void get_build_id(char* build_id)
{
   // Dummy string is replaced by script "generate_userblock.sh"
   strcpy(build_id,"DUMMY_INSERT_BUILD_ID_HERE");
}


/* Copies contents of file "infilename" to "outfilename". */
void copy_userblock(char* outfilename, char* infilename)
{
   FILE* fout = fopen(outfilename,"rb+");
   FILE* fin  = fopen(infilename, "r");
   rewind(fout); // Set file pointer to first line
   int c;
   // Copy line for line until EOF is encountered
   do {
      c = fgetc (fin);
      if (c != EOF) fputc((char) c, fout);
   } while (c != EOF);
   fclose(fin);
   fclose(fout);
}


/* Creates and dumps the userblock into the file "filename".
 *
 * This routine creates the userblock and dumps it into the file
 * "filename". For this, it first copies the contents of the parameter
 * file "inifilename" into the output file "filename". After this, the
 * code information (Git Hash etc.) is appended and, finally, the build
 * information (CMake configuration etc.). The latter two are generated
 * and inserted here during the compile process using the
 * "generate_userblock.sh" script.
 */
void print_userblock(char* filename, char* inifilename)
{
   FILE* fp = fopen(filename, "w");
   rewind(fp);

   fprintf(fp, "{[( START USERBLOCK )]}\n");

   // Dump parameter file section into userblock
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

   /* the shell script "generate_userblock.sh" copies
      this file into the build directory and inserts
      print commands for the code and build info here. */

   //INSERT_CODE_INFO_HERE

   //INSERT_BUILD_INFO_HERE

   fprintf(fp, "{[( END USERBLOCK )]}\n");

   fclose(fp);
}


/* Writes code ID and code information to "filename". Both are
 * determined during compile time by the "generate_userblock.sh"
 * script.
 */
void print_code_info(char* filename)
{

   FILE* fp = fopen(filename, "w");
   rewind(fp);

   fprintf(fp, "{[( CODE ID )]}\n");
   // Dummy string is replaced by script "generate_userblock.sh"
   fprintf(fp, "DUMMY_INSERT_CODE_ID_HERE\n");

   /* the shell script generate_userblock.sh copies
      this file into the build directory and inserts
      print commands for the code info here. */

   //INSERT_CODE_INFO_HERE

   fclose(fp);
}


/* Writes build and code ID as well as code information to "filename".
 * Both are determined during compile time by the
 * "generate_userblock.sh" script.
 */
void print_build_info(char* filename)
{

   FILE* fp = fopen(filename, "w");
   rewind(fp);

   fprintf(fp, "{[( BUILD ID )]}\n");
   // Dummy string is replaced by script "generate_userblock.sh"
   fprintf(fp, "DUMMY_INSERT_BUILD_ID_HERE\n");

   fprintf(fp, "{[( CODE ID )]}\n");
   // Dummy string is replaced by script "generate_userblock.sh"
   fprintf(fp, "DUMMY_INSERT_CODE_ID_HERE\n");

   /* the shell script "generate_userblock.sh" copies
      this file into the build directory and inserts
      print commands for the build info here. */

   //INSERT_BUILD_INFO_HERE

   fclose(fp);
}
