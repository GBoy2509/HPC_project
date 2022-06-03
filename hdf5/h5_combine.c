/* 
 *  This example illustrates how to write and read data in an existing
 *  dataset.  It is used in the HDF5 Tutorial.
 */

#include "stdio.h"
#include "stdlib.h"
#include "hdf5.h"
#define FILE "combine.h5"

int main() {

  hid_t       file_id, dataset_id, group_id, dataspace_id;  /* identifiers */
  hsize_t     dims[3];
  herr_t      status;
  int         i, j, k,  dset_data[5][6][4];
  int * vec1 = (int*)malloc(5*6*4*sizeof(int));   
          free(vec1);

    for (j = 0; j < 5; j++) {
    for (i = 0; i < 6; i++){
      for (k = 0; k < 4; k++)
        dset_data[j][i][k] = 200 + 100*j + 10*i + k;
    }
  }  

  /* Create a new file using default properties. */
  file_id = H5Fcreate("combine.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  dims[0] = 5; 
  dims[1] = 6;
  dims[2] = 4;
  dataspace_id = H5Screate_simple(3, dims, NULL);

  /* Create a group named "/MyGroup" in the file. */
  group_id = H5Gcreate2(file_id, "/MyGroup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Create the dataset. */
  dataset_id = H5Dcreate2(file_id, "/MyGroup/IntArray", H5T_STD_I32BE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


  //  printf("original dset_data[0][0][0]:%2d\n", dset_data[0][0][0]);

   /* Write the first dataset. */
  status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     dset_data);

   /* Close the data space for the first dataset. */
  status = H5Sclose(dataspace_id);

   /* Close the first dataset. */
  status = H5Dclose(dataset_id);
   /* Close the group. */
  status = H5Gclose(group_id);
  status = H5Fclose(file_id);
}
