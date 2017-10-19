# each analyze-function takes the output-lines of a flexi-run 
import math

# extract the L2 error of the last timestep
def get_last_L2_error(lines) :
   for l in lines[-15:] :
      if "L_2" in l :
         tmp = l.split(":")[1]
   return [float(x) for x in tmp.split()]

# extract the L_inf error of the last timestep
def get_last_Linf_error(lines) :
   for l in lines[-15:] :
      if "L_inf" in l :
         tmp = l.split(":")[1]
         return [float(x) for x in tmp.split()]

def get_last_number(lines) :
   for line in reversed(lines) :
      tmp = line.split(' ')
      for t in reversed(tmp) :
         try :
            return float(t)
         except :
            pass

def get_cpu_per_dof(lines) :
   for line in reversed(lines) :
      if "CALCULATION TIME PER TSTEP/DOF: [" in line :
         return float(line.split("[")[1].split("sec")[0])

def calcOrder_h(h,E) :
    if len(h) != len(E) :
        return -1

    order = []
    for i in range(1,len(h)) :
        dh=1.0/(h[i]/h[i-1])
        dE=E[i]/E[i-1]
        order.append(math.log(dE)/math.log(dh))

    return order
