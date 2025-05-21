#include <fenv.h>

void enableFPE() {
    feclearexcept(FE_ALL_EXCEPT);
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}

void disableFPE() {
    feclearexcept(FE_ALL_EXCEPT);
    fedisableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}
