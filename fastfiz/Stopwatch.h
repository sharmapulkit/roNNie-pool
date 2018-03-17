/*
   This class determines how much time has elapsed since the last time
   reset was called.
   */

#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_

#include <time.h>

typedef struct timeval timeval;

/** Utility class heirarchy to measure time in various ways.
 *  Base class implements general time handling with arbitrary units.
 *  The virtual now() method should be implemented in subclasses.
 */
class Stopwatch {
public:
  /** Create a new stopwatch. Not running. */
  Stopwatch() : _extra(0.0), _running(false){};
  /** Reset stopwatch to zero and stop. */
  void clear() {
    _extra = 0.0;
    _running = false;
  };
  /** Reset stopwatch to zero and start. */
  void restart() {
    _extra = 0.0;
    _start = now();
    _running = true;
  };
  /** Get current reading on stopwatch. */
  double getElapsed() const {
    return (_running ? now() - _start : 0) + _extra;
  };
  /** Stop the stopwatch. Retain time reading. */
  void stop() {
    if (_running)
      _extra += now() - _start;
    _running = false;
  };
  /** Start the stopwatch. Retain time reading. */
  void start() {
    if (!_running) {
      _start = now();
      _running = true;
    };
  };
  /** Create a new stopwatch of same type (factory). */
  virtual Stopwatch *create() const = 0;

protected:
  /** Get the current time in stopwatch units */
  virtual double now() const = 0;

private:
  /** Time when stopwatch is counting from (if running) */
  double _start;
  /** Extra time to add (previous laps) */
  double _extra;
  /** Is the stopwatch running */
  bool _running;
};

/** A stopwatch that measures in actual seconds, but not necessarily wall time
 */
class TimevalStopwatch : public Stopwatch {
protected:
  /** Get current time in seconds since the epoch */
  virtual double now() const;
  /** Get current time as a standard timeval */
  virtual timeval now_timeval() const = 0;
};

/** A stopwatch that measures wall time */
class RealTimeStopwatch : public TimevalStopwatch {
public:
  virtual Stopwatch *create() const { return new RealTimeStopwatch(); };

protected:
  virtual timeval now_timeval() const;
};

/** A stopwatch that measures in CPU time (ignoring time spent by system or
 * other processes) */

class CPUStopwatch : public TimevalStopwatch {
public:
  virtual Stopwatch *create() const { return new CPUStopwatch(); };

protected:
  virtual timeval now_timeval() const;
};

/** A stopwatch measuring virtual time kept in a static counter.
 *  Each tick is 3e-3 time units.
 */
class VirtualStopwatch : public Stopwatch {
public:
  virtual Stopwatch *create() const { return new VirtualStopwatch(); };
  static void step() { counter++; };

protected:
  virtual double now() const { return counter * 3e-3; };

private:
  static unsigned long counter;
};

/** A stopwatch measuring virtual time kept in a static counter.
 *  Each tick is 1e-4 time units.
 */
class VirtualStopwatch2 : public Stopwatch {
public:
  virtual Stopwatch *create() const { return new VirtualStopwatch2(); };
  static void step() { counter++; };

protected:
  virtual double now() const { return counter * 1e-4; };

private:
  static unsigned long counter;
};

#endif
