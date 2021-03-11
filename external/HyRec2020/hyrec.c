/*************************************************************************/
/*                    HYREC-2 MAIN FUNCTION                              */
/*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "history.h"

int main(void) {
  remove("output_xe.dat");
  
  HYREC_DATA rec_data;

  double zmax = 8000.;
  double zmin = 0.;

  rec_data.path_to_hyrec = "";
  hyrec_allocate(&rec_data, zmax, zmin);
  rec_get_cosmoparam(stdin, stderr, rec_data.cosmo);
  hyrec_compute(&rec_data, MODEL);
  if (rec_data.error == 1) fprintf(stderr,"%s\n",rec_data.error_message);
  else {
    double z = zmax;
    char file[200];
    FILE *fout;
    fout = fopen("output_xe.dat", "a");
    while (z > zmin) {
      fprintf(fout, "%f %1.10E %1.10E\n",z,hyrec_xe(z, &rec_data),hyrec_Tm(z, &rec_data));
      z -= 1.;
      }
    fclose(fout);
  }
  hyrec_free(&rec_data);
  return 0;

}
