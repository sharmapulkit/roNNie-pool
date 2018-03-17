#include "Stopwatch.h"
#include <sys/resource.h>
#include <sys/time.h>

double TimevalStopwatch::now() const {
  timeval tv = now_timeval();
  return tv.tv_sec + tv.tv_usec / (double)1e6;
}

timeval RealTimeStopwatch::now_timeval() const {
  timeval tv;
  gettimeofday(&tv, NULL);
  return tv;
}

timeval CPUStopwatch::now_timeval() const {
  rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return ru.ru_utime;
}

unsigned long VirtualStopwatch::counter = 0;
unsigned long VirtualStopwatch2::counter = 0;
