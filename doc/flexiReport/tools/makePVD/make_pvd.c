#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int main(int argc, char *argv[]){

	if(argc != 2){

		fprintf(stderr, "Usage: %s <file_increment>\n", argv[0]);
		exit(1);

	}

	bool exvar = true;
	FILE *test_file;
	FILE *output = fopen("einfeldt.pvd", "w");
    double nvar = 0;
	char filename[1000];

	double incr = atof(argv[1]);

	fprintf(output, "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
    fprintf(output, "  <Collection>\n");

	
    

	do {

		sprintf(filename, "data/leblanc_Solution_%017.9f.vtm", nvar);

		test_file = fopen(filename, "r");

		if(!test_file) exvar = false;
		else {
			fclose(test_file);
		    printf("Writing file: %s\n", filename);
			/* write it into the pvd file */
            fprintf(output, "     <DataSet timestep=\"%f\" part=\"0\" file=\"%s\"/>\n", nvar, filename);
		    nvar += incr;
		}

	} while (exvar);

    fprintf(output, "  </Collection>\n</VTKFile>\n");


}
