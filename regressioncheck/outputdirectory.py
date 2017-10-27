import os

class OutputDirectory() :
    output_dir = "output_dir"
    def __init__(self, parent, name, number = -1, mkdir=True) :

        self.number = number
        self.parent = parent

        # set parent directory for subfolder creation
        if self.parent :
            parent_dir = self.parent.target_directory
        else :
            parent_dir = OutputDirectory.output_dir 

        # numbering of directory (if a number is supplied)
        if number >= 0 :
            self.target_directory = os.path.join(parent_dir, "%s_%04d" %(name, number))
        else :
            self.target_directory = os.path.join(parent_dir, name)

        # create directory if it is non-existent
        if mkdir :
            if not os.path.exists(self.target_directory) :
                os.makedirs(self.target_directory)


