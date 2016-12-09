/*
!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
*/

#include "visu3DReader.h"
#include "../posti_plugin.h"

#include "hdf5.h"
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkHexahedron.h>
#include <vtkQuad.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include "vtkMultiBlockDataSet.h"
#include <vtkUnstructuredGrid.h>


#include <libgen.h>
#include <unistd.h>
#include <algorithm>

// MPI
#include "vtkMultiProcessController.h"
vtkStandardNewMacro(visu3DReader);
vtkCxxSetObjectMacro(visu3DReader, Controller, vtkMultiProcessController);

/*
 * Construtor of State Reader
 */
visu3DReader::visu3DReader()
{
   this->FileName = NULL;
   this->InputNsuper = 0;
   this->NodeTypeVisu = NULL;
   this->Mode2d = 0;
   this->ParameterFileOverwrite = NULL;
   this->MeshFileOverwrite = NULL;
   this->SetNumberOfInputPorts(0);


   // Setup the selection callback to modify this object when an array
   // selection is changed.
   // Used to tell the visu3DReader, that we (un)selected a state,primite or derived quantity
   // and that the 'Apply' button becomes clickable to reload the data (load the selected quantities)

   this->SelectionObserver = vtkCallbackCommand::New();
   this->SelectionObserver->SetCallback(&visu3DReader::SelectionModifiedCallback);
   this->SelectionObserver->SetClientData(this);
   // create array for state,primitive and derived quantities
   this->VarDataArraySelection = vtkDataArraySelection::New();
   // add an observer 
   this->VarDataArraySelection->AddObserver(vtkCommand::ModifiedEvent, this->SelectionObserver);
}

/*
 * This function is called when a file is inserted into the Pipeline-browser.
 * Here we load the variable names (state,primitive, derived).
 */
int visu3DReader::RequestInformation(vtkInformation *,
      vtkInformationVector **,
      vtkInformationVector *outputVector)
{
   vtkInformation* info = outputVector->GetInformationObject(0);
   info->Set(vtkAlgorithm::CAN_HANDLE_PIECE_REQUEST(), 1);
   //#if USE_MPI
   // Set up MPI communicator   
   this->Controller = NULL;
   this->SetController(vtkMultiProcessController::GetGlobalController());
   if (this->Controller == NULL) {
      NumProcesses = 1;
      ProcessId = 0;
   } else {
      NumProcesses = this->Controller->GetNumberOfProcesses();
      ProcessId = this->Controller->GetLocalProcessId();
   }

   vtkMPICommunicator *communicator = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
   mpiComm = MPI_COMM_NULL;
   if (communicator) {
      mpiComm = *(communicator->GetMPIComm()->GetHandle());
   }

   // get the info object of the output ports 
   vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);

   // sets the number of pieces to the number of processsors
   outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

   // if we have more then one file loaded at once (timeseries)
   // we have to set the number and range of the timesteps
   double timeRange[2] = {Timesteps[0], Timesteps[Timesteps.size()-1]};
   outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &Timesteps[0], Timesteps.size());
   outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

   // We take the first state file and use it to read the varnames
   if(ProcessId == 0) std::cout << "RequestInformation: State file: " << FileNames[0] << std::endl;

   int fcomm;
   fcomm = MPI_Comm_c2f(mpiComm);
   MPI_Barrier(mpiComm);
   struct CharARRAY varnames;
   // Call Posti-function requestInformation:
   // This function returns the varnames of state, primitive and derived quantities
   int str_len = strlen(FileNames[0].c_str());
   __mod_visu3d_MOD_visu3d_requestinformation(&fcomm, &str_len, FileNames[0].c_str(), &varnames);

   MPI_Barrier(mpiComm);
   // We copy the varnames to the corresponding DataArraySelection objects.
   // These objects are used to build the gui. 
   // (see the functions below: 
   //   DisableAllVar...Arrays, EnableAllVar...Arrays, GetNumberOfVar...Arrays, GetVar...ArrayName, GetVar...ArrayStatus, SetVar...ArrayStatus,
   // )
   for (int iVar=0; iVar<varnames.len/255; iVar++) {
      char tmps[255];
      strncpy(tmps, varnames.data+iVar*255, 255);
      std::string varname(tmps);
      varname = varname.substr(0,varname.find(" "));

      if (!this->VarDataArraySelection->ArrayExists(varname.c_str())) {
         // function DisableArray deselects the varname in the gui and if not existend inserts the varname         
         this->VarDataArraySelection->DisableArray(varname.c_str());


         // Select Density, FV_Elems by default
         if (varname.compare("Density") == 0) this->VarDataArraySelection->EnableArray(varname.c_str());
         if (varname.compare("FV_Elems") == 0) this->VarDataArraySelection->EnableArray(varname.c_str());
      }
   }
   return 1; 
}

/*
 * This function is called whenever a filename is loaded.
 * If multiple files are selected in the file-dialog, this function is called multiple times.
 * Attention: For multiple files, we assume a timeseries.
 */
void visu3DReader::AddFileName(const char* filename_in) {
   // append the filename to the list of filenames
   this->FileNames.push_back(filename_in);
   this->Modified();

   // open the file with HDF5 and read the attribute 'time' to build a timeseries  
   hid_t state = H5Fopen(filename_in, H5F_ACC_RDONLY, H5P_DEFAULT);
   hid_t attr = H5Aopen(state, "Time", H5P_DEFAULT);
   if(ProcessId == 0) std::cout<<"attribute Time"<<attr<<"\n";
   double time; 
   if (attr > -1){
      hid_t attr_type = H5Aget_type( attr );
      H5Aread(attr, attr_type, &time);
      Timesteps.push_back(time);
   }else{
      Timesteps.push_back(0.);
   }
   H5Aclose(attr);
   H5Fclose(state);
}

void visu3DReader::RemoveAllFileNames() {
   this->FileNames.clear();
   this->Timesteps.clear();
}

// Get number of state file in the filenames list closest to the requested time
int visu3DReader::FindClosestTimeStep(double requestedTimeValue)
{
   int ts = 0;
   double mindist = fabs(Timesteps[0] - requestedTimeValue);

   for (unsigned int i = 0; i < Timesteps.size(); i++) {
      double tsv = Timesteps[i];
      double dist = fabs(tsv - requestedTimeValue);
      if (dist < mindist) {
         mindist = dist;
         ts = i;
      }
   }
   return ts;
}

//void visu3DReader::SetNodeTypeVisu(const char* nodetypevisu) {
   //NodeTypeVisu = nodetypevisu;
//}

vtkStringArray* visu3DReader::GetNodeTypeVisuList() {
   vtkStringArray* arr = vtkStringArray::New();
   arr->InsertNextValue("VISU");
   arr->InsertNextValue("GAUSS");
   arr->InsertNextValue("GAUSS-LOBATTO");
   arr->InsertNextValue("VISU_INNER");
   return arr;
};

/*
 * This function is called, when the user presses the 'Apply' button.
 * Here we call the Posti and load all the data.
 */
   int visu3DReader::RequestData(
         vtkInformation *vtkNotUsed(request),
         vtkInformationVector **vtkNotUsed(inputVector),
         vtkInformationVector *outputVector)
{

   // get the number of the timestep to load
   int timestepToLoad = 0;
   vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
   if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
      // get the requested time
      double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
      timestepToLoad = FindClosestTimeStep(requestedTimeValue);

      // get the output and set the time value
      vtkSmartPointer<vtkUnstructuredGrid> output = vtkUnstructuredGrid::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
      //output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), Timesteps[timestepToLoad]);
   }

   vtkDataObject* doOutput = outInfo->Get(vtkDataObject::DATA_OBJECT());
   vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::SafeDownCast(doOutput);
   if (!mb)
   {
      std::cout << "DownCast to MultiBlockDataset Failed!" << std::endl;
      return 0;
   }


   // Change to directory of state file (path of mesh file is stored relative to path of state file)
   char dir[255];
   strcpy(dir, FileNames[timestepToLoad].c_str());
   int ret = chdir(dirname(dir));
   if (ret != 0 ) {
      if(ProcessId == 0) std::cout << "dirfolder of statefile not found: " << dirname(dir) << " \n"; 
   }

   // Write temporary parameter file for Posti tool 
   // TODO: only root-proc
   //if (ProcessId == 0) {

   // get temporary file for the posti parameter files
   setlocale(LC_ALL, "C"); 
   char posti_filename[255];
   char prm_filename[255];
   std::strcpy(posti_filename, "/tmp/f2p_posti_XXXXXX.ini");
   int posti_unit = mkstemps(posti_filename,4);


   if (strlen(ParameterFileOverwrite) > 0) {
      // use parameter file to overwrite userblock-parameterfile 
      std::strcpy(prm_filename, ParameterFileOverwrite);
   } else {
      std::strcpy(prm_filename, "");
   }

   dprintf(posti_unit, "NVisu = %d\n", InputNsuper); // insert NVisu
   dprintf(posti_unit, "NodeTypeVisu = %s\n", NodeTypeVisu); // insert NodeType
   dprintf(posti_unit, "VisuDimension = %s\n", (this->Mode2d ? "2" : "3"));
   if (strlen(MeshFileOverwrite) > 0) {
      dprintf(posti_unit, "MeshFile = %s\n", MeshFileOverwrite);
   }

   int totalVars = 0;
   // write selected state varnames to the parameter file
   int nVars = VarDataArraySelection->GetNumberOfArrays();
   for (int i = 0; i< nVars; ++i)
   {
      const char* name = VarDataArraySelection->GetArrayName(i);
      if (VarDataArraySelection->ArrayIsEnabled(name)) {
         dprintf(posti_unit, "VarName = %s\n", name) ;
         totalVars++;
      }
   }
   close(posti_unit);

   int fcomm = MPI_Comm_c2f(mpiComm);
   MPI_Barrier(mpiComm);
   // call Posti tool 
   // the arrays coords_*, values_* and nodeids_* are allocated in the Posti tool 
   // and contain the vtk data
   int strlen_prm = strlen(prm_filename);
   int strlen_posti = strlen(posti_filename);
   int strlen_state = strlen(FileNames[timestepToLoad].c_str()); 
   __mod_visu3d_MOD_visu3d_cwrapper(&fcomm, 
         &strlen_prm, prm_filename, 
         &strlen_posti, posti_filename, 
         &strlen_state, FileNames[timestepToLoad].c_str(),
         &coords_DG,&values_DG,&nodeids_DG,
         &coords_FV,&values_FV,&nodeids_FV,&varnames,&components);
   MPI_Barrier(mpiComm);

   this->Blocks.resize(0);

   // insert data into output (DG)
   this->Blocks.resize(this->Blocks.size()+1);
   vtkSmartPointer<vtkUnstructuredGrid> output_DG = this->Blocks[this->Blocks.size()-1];
   if (!output_DG)
   {
      output_DG = vtkUnstructuredGrid::New();
      this->Blocks[this->Blocks.size()-1] = output_DG;
      output_DG->Delete();
   }
   // Insert data into output
   InsertData(output_DG, &coords_DG, &values_DG, &nodeids_DG, &varnames, &components);

   // insert data into output (FV)
   this->Blocks.resize(this->Blocks.size()+1);
   vtkSmartPointer<vtkUnstructuredGrid> output_FV = this->Blocks[this->Blocks.size()-1];
   if (!output_FV)
   {
      output_FV = vtkUnstructuredGrid::New();
      this->Blocks[this->Blocks.size()-1] = output_FV;
      output_FV->Delete();
   }
   // Insert data into output
   InsertData(output_FV, &coords_FV, &values_FV, &nodeids_FV, &varnames, &components);


   // the nodeids are copied in the 'InsertData' function. Therefore deallocate the arrays in the Posti tool
   __mod_visu3d_MOD_visu3d_dealloc_nodeids();
   MPI_Barrier(mpiComm);
   int minusBlock = 0;


   if(ProcessId == 0) std::cout << "minusBlock " << minusBlock << " " << std::endl;
   int firstBlock=0; 
   mb->SetNumberOfBlocks(this->Blocks.size()-minusBlock); // missing -1
   for(unsigned int i=0; i<this->Blocks.size()-minusBlock; i++) // missing -1
   {
      vtkUnstructuredGrid* nthOutput = this->Blocks[i];
      firstBlock=firstBlock+1;
      mb->SetBlock(i, nthOutput);
   }


   this -> Modified(); // tell paraview to render data
   MPI_Barrier(mpiComm);

   return 1;
}


/*
 * This function inserts the data, loaded by the Posti tool, into a ouput
 */
void visu3DReader::InsertData(vtkSmartPointer<vtkUnstructuredGrid> &output, struct DoubleARRAY* coords,
      struct DoubleARRAY* values, struct IntARRAY* nodeids, struct CharARRAY* varnames, struct IntARRAY* components) {

   // create points(array), cellarray and a hex;
   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
   vtkSmartPointer<vtkCellArray> cellarray = vtkSmartPointer<vtkCellArray>::New();
   // create a 3D double array (must be 3D even we use a 2D Posti tool, since paraview holds the data in 3D)
   vtkSmartPointer <vtkDoubleArray> pdata = vtkSmartPointer<vtkDoubleArray>::New();
   pdata->SetNumberOfComponents(3); // 3D

   // assign the coords-data (loaded by the Posti tool) to the general pdata array 
   // (no copy of data!!! pdata uses the coords-data array)

   pdata->SetArray(coords->data, coords->len, 1);
   // now we use this vtkDouble Array and assign it to the points-array
   // (again no copy of data!!!)
   points->SetData(pdata);
   points->SetNumberOfPoints(coords->len/3);
   // set the points array to be used for the output
   output->SetPoints(points);


   if (this->Mode2d) {
      vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
      // Use the nodeids to build quads
      // (here we must copy the nodeids, we can not just assign the array of nodeids to some vtk-structure)
      int gi = 0;
      // loop over all Quads
      for (int iQuad=0; iQuad<nodeids->len/4; iQuad++) {
         // each Quad has 4 points
         for (int i=0; i<4; i++) {
            quad->GetPointIds()->SetId(i, nodeids->data[gi]);
            gi++;
         }
         // insert the quad into the cellarray
         cellarray->InsertNextCell(quad);
      }
      output->SetCells({VTK_QUAD}, cellarray);
   } else {
      vtkSmartPointer<vtkHexahedron> hex = vtkSmartPointer<vtkHexahedron>::New();
      // Use the nodeids to build hexas
      // (here we must copy the nodeids, we can not just assign the array of nodeids to some vtk-structure)
      int gi = 0;
      // loop over all Hexas
      for (int iHex=0; iHex<nodeids->len/8; iHex++) {
         // each Hex has 8 points
         for (int i=0; i<8; i++) {
            hex->GetPointIds()->SetId(i, nodeids->data[gi]);
            gi++;
         }
         // insert the hex into the cellarray
         cellarray->InsertNextCell(hex);
      }
      output->SetCells({VTK_HEXAHEDRON}, cellarray);
   }

   // assign the actual data, loaded by the Posti tool, to the output
   std::vector<vtkSmartPointer<vtkDoubleArray> > dataArrays;
   unsigned int nVarCombine = varnames->len/255;
   if (nVarCombine > 0) {
      unsigned int nVar = 0;
      for (unsigned int iVar=0;iVar<nVarCombine;iVar++) {
         nVar += components->data[iVar];
      }
      unsigned int sizePerVar = values->len/nVar;
      int dataPos = 0;
      // loop over all loaded variables
      for (unsigned int iVar = 0; iVar < nVarCombine; iVar++) {
         // For each variable, create a new array and set the number of components to 1
         // Each variable is loaded separately. 
         // One might implement vector variables (velocity), but then must set number of componenets to 2/3
         dataArrays.push_back( vtkSmartPointer<vtkDoubleArray>::New() );
         dataArrays.at(iVar)->SetNumberOfComponents(components->data[iVar]);
         // Assign data of the variable to the array.
         // (no copy of data!!!)
         dataArrays.at(iVar)->SetArray(values->data+dataPos, sizePerVar*components->data[iVar], 1);
         dataPos += sizePerVar * components->data[iVar];
         // set name of variable
         char tmps[255];
         strncpy(tmps, varnames->data+iVar*255, 255);
         std::string varname(tmps);
         varname = varname.substr(0,varname.find(" "));
         dataArrays.at(iVar)->SetName(varname.c_str());
         // insert array of variable into the output
         output->GetPointData()->AddArray(dataArrays.at(iVar));
      }
   }
}

visu3DReader::~visu3DReader(){
   delete [] FileName;

   // Finalize the Posti tool (deallocate all arrays)
   __mod_visu3d_MOD_visu3d_finalize();

   this->VarDataArraySelection->Delete();
}

/*
 * The following functions create the interaction between this Reader and 
 * the gui, which is defined in the visu2DReader.xml
 * They return the number of available variables (state,primitive, derived)
 * and return the names of the variables, ....
 */

void visu3DReader::DisableAllVarArrays()
{
   this->VarDataArraySelection->DisableAllArrays();
}
void visu3DReader::EnableAllVarArrays()
{
   this->VarDataArraySelection->EnableAllArrays();
}
int visu3DReader::GetNumberOfVarArrays()
{
   return this->VarDataArraySelection->GetNumberOfArrays();
}

const char* visu3DReader::GetVarArrayName(int index)
{
   if (index >= ( int ) this->GetNumberOfVarArrays() || index < 0)
   {
      return NULL;
   }
   else
   {
      return this->VarDataArraySelection->GetArrayName(index);
   }
}
int visu3DReader::GetVarArrayStatus(const char* name)
{
   return this->VarDataArraySelection->ArrayIsEnabled(name);
}

void visu3DReader::SetVarArrayStatus(const char* name, int status)
{
   if (status)
   {
      this->VarDataArraySelection->EnableArray(name);
   }
   else
   {
      this->VarDataArraySelection->DisableArray(name);
   }
}

void visu3DReader::SelectionModifiedCallback(vtkObject*, unsigned long,
      void* clientdata, void*)
{
   static_cast<visu3DReader*>(clientdata)->Modified();
}

int visu3DReader::FillOutputPortInformation(
      int vtkNotUsed(port), vtkInformation* info)
{
   info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
   return 1;
}


