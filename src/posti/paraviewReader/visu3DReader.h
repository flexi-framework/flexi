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

#ifndef VISU3DREADER_H
#define VISU3DREADER_H

#include <vtkUnstructuredGridAlgorithm.h>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <vector>

#include "vtkMPI.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIController.h"

#include "vtkDataArraySelection.h"
#include "vtkCallbackCommand.h"
#include <vtkSmartPointer.h>

#include "../posti_pluginTypes.h"

#include "vtkIOParallelModule.h" // For export macro
#include "vtkMultiBlockDataSetAlgorithm.h"

// MPI
class vtkMultiProcessController;
// MPI

class VTKIOPARALLEL_EXPORT visu3DReader :  public vtkMultiBlockDataSetAlgorithm
{
   public:
      vtkTypeMacro(visu3DReader,vtkMultiBlockDataSetAlgorithm);

      static visu3DReader *New();

      // macros to set Filename, InputNsuper (see visu2DReader.xml)
      // gui interaction
      vtkSetStringMacro(FileName);
      vtkSetMacro(InputNsuper,int);
      vtkSetMacro(Mode2d,int);

      // Adds names of files to be read. The files are read in the order they are added.
      void AddFileName(const char* fname);

      // Remove all file names.
      void RemoveAllFileNames();

      vtkSetStringMacro(ParameterFileOverwrite);

      int GetNumberOfVarArrays();
      const char* GetVarArrayName(int index);
      int GetVarArrayStatus(const char* name);
      void SetVarArrayStatus(const char* name, int status);
      void DisableAllVarArrays();
      void EnableAllVarArrays();

      // MPI
      void SetController(vtkMultiProcessController *);
      vtkGetObjectMacro(Controller, vtkMultiProcessController);
      // MPI

      void ConvertToFortran(char* fstring, const char* cstring);

      // struct to exchange arrays between fortran and C      
      struct DoubleARRAY coords_DG;
      struct DoubleARRAY values_DG;
      struct IntARRAY  nodeids_DG;
      struct DoubleARRAY coords_FV;
      struct DoubleARRAY values_FV;
      struct IntARRAY  nodeids_FV;
      struct CharARRAY varnames;
      struct IntARRAY components;

      int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
      int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

      void InsertData(vtkSmartPointer<vtkUnstructuredGrid> &output, struct DoubleARRAY* coords,
            struct DoubleARRAY* values, struct IntARRAY* nodeids, struct CharARRAY* varnames, struct IntARRAY* components);

      vtkDataArraySelection* VarDataArraySelection;

      // The observer to modify this object when the array selections are modified.
      vtkCallbackCommand* SelectionObserver;

      // Callback registered with the SelectionObserver.
      static void SelectionModifiedCallback(vtkObject* caller, unsigned long eid,
            void* clientdata, void* calldata);


      char* FileName;
      int InputNsuper;   // NVisu (see visu2DReader.xml)
      int Mode2d;
      char* ParameterFileOverwrite;

      int NumProcesses;
      int ProcessId;

      // all loaded filenames, timesteps (multiple for timeseries)
      std::vector<std::string> FileNames;
      std::vector<double> Timesteps;

      int FindClosestTimeStep(double requestedTimeValue);

   protected:
      visu3DReader();
      ~visu3DReader();

      virtual int FillOutputPortInformation(int port, vtkInformation* info);

      std::vector<vtkSmartPointer<vtkUnstructuredGrid> > Blocks;


   private:
      // MPI
      vtkMultiProcessController *Controller;
      MPI_Comm mpiComm;
};

#endif //VISU3DREADER_H
