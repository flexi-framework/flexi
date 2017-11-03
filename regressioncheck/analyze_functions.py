import math

def get_last_L2_error(lines) :
   """Get L_2 eror value from a set of lines for the last timestep.
   The set of lines correspond to the output-lines of a flexi-run"""
   for l in lines[-15:] :
      if "L_2" in l :
         tmp = l.split(":")[1]
   return [float(x) for x in tmp.split()]

def get_last_Linf_error(lines) :
   """Get L_inf eror value from a set of lines for the last timestep
   The set of lines correspond to the output-lines of a flexi-run"""
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
   """Get the PID value from a set of lines
   The set of lines correspond to the output-lines of a flexi-run"""
   for line in reversed(lines) :
      if "CALCULATION TIME PER TSTEP/DOF: [" in line :
         return float(line.split("[")[1].split("sec")[0])

def calcOrder_h(h,E) :
    """Determine the order of convergence for a list of grid spacings h and errors E"""
    h = [float(elem) for elem in h]
    E = [float(elem) for elem in E]
    if len(h) != len(E) :
        return -1

    order = []
    for i in range(1,len(h)) :
        dh=1.0/(h[i]/h[i-1])
        if E[i-1] == 0.0 :
            order.append(0.0)
        else :
            dE=E[i]/E[i-1]
            order.append(math.log(dE)/math.log(dh))

    return order

def calcOrder_p(p,E) :
    """Determine the order of convergence for a list of polynomial degrees p and errors E"""
    p = [float(elem) for elem in p]
    E = [float(elem) for elem in E]
    if len(p) != len(E) :
        return -1

    order = []
    for i in range(1,len(p)) :
        dp=1.0/((p[i]+1.0)/(p[i-1]+1.0))
        if E[i-1] == 0.0 :
            order.append(0.0)
        else :
            dE=E[i]/E[i-1]
            order.append(math.log(dE)/math.log(dp))

    return order
