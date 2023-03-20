#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "astroFns.h"
#include "structures.h"

#include "info.h"
#include "mathFns.h"
#include "random.h"

void applyWeather(struct obsfilekeywords World[], struct filekeywords *Paramfile, long *idum);
int getWeather(double PrClearnight, long *idum);
void loadWeatherProfile(char profileName[], double WeatherProfile[], struct filekeywords *Paramfile);
