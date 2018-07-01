cdef extern from "nem_exe.h":
   cpdef int nem(const char* Fname,
        const int nk,
        const char* algo,
        const float beta,
        const char* convergence,
        const float convergence_th,
        const char* format,
        const int it_max,
        const int dolog,
        const char* model_family,
        const char* proportion,
        const char* dispersion,
        const int init_mode);
