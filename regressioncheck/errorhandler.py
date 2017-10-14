import os

class ErrorHandler() :
    successful = True
    errors = []

    def __init__(self, parent, name, number = -1, mkdir=True) :
        self.number = number
        self.parent = parent
        if self.parent :
            parent_dir = self.parent.directory
        else :
            parent_dir = "reggie_outdir"
        if number >= 0 :
            self.directory = os.path.join(parent_dir, "%s_%04d" %(name, number))
        else :
            self.directory = os.path.join(parent_dir, name)

        if mkdir :
            if not os.path.exists(self.directory) :
                os.mkdir(self.directory)  # create example directory
