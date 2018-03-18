#ifndef _LogFile_h
#define _LogFile_h

#ifdef SWIG
%include exception.i
%include "std_string.i"
%{
#include "LogFile.h"
%}
#endif /* SWIG */

#include "FastFiz.h"
#include "Noise.h"
#include "Rules.h"
#include <fstream>
#include <iostream>
#include <string>

namespace Pool {

/**
 * Output game/table states and shots to a stream in a text-based format
 * supported by web-based game viewer.
 *
 * This class allows for testing shots offline without writing a complete AI.
 * The resulting log file can
 * then be used to visualize the results of the shots.
 *
 * Log files can also be used as a convenient way to store shots or table states
 * for later use.
 */
class LogWriter {
public:
  /**
   * Start a new log.
   *
   * \param os Output stream to write log.
   * \param gameType Type of game being played (e.g. GT_EIGHTBALL).
   * \param noise Noise object being used. Noise info will be written to log.
   * \param agentName Name of agent executing the shots (usually ai.getName()).
   * \param opponentName Name of opponent to be written to log.
   */
  LogWriter(std::ostream &os, const GameType gameType = GT_NONE,
            const Noise *const noise = 0,
            const std::string agentName = "Default Agent",
            const std::string opponentName = "")
      : _agentName(agentName), _opponentName(opponentName), _noise(noise),
        _ofs(NULL), _os(os) {
    writeHeader(gameType);
  };

  /**
   * Start a new log file (conveinent version of constructor for files).
   *
   * \param filename Name of file to write.
   * \param gameType Type of game being played (e.g. GT_EIGHTBALL).
   * \param noise Noise object being used. Noise info will be written to log.
   * \param agentName Name of agent executing the shots (usually ai.getName()).
   * \param opponentName Name of opponent to be written to log.
   */
  LogWriter(const char *const filename, const GameType gameType = GT_NONE,
            const Noise *const noise = 0,
            const std::string agentName = "Default Agent",
            const std::string opponentName = "")
      : _agentName(agentName), _opponentName(opponentName), _noise(noise),
        _ofs(new std::ofstream(filename)), _os(*_ofs) {
    writeHeader(gameType);
  };

  /**
   * Set Noise type for future shots. Useful when noise is variable.
   */
  void setNoise(const Noise *const noise) { _noise = noise; };

  /**
   * Set agent names.
   * \param agentName Name of agent executing future shots.
   * \param opponentName Name of the other agent.
   */
  void setAgents(const std::string agentName, const std::string opponentName = "") {
    _agentName = agentName;
    _opponentName = opponentName;
  };

  /**
   * Make opponent agent active (for future shots) and set active agent as
   * opponent.
   */
  void swapAgents() { setAgents(_opponentName, _agentName); };

  /**
   * Add Entry to log file with given game state.
   */
  void write(const GameState &gs);

  /**
   * Add Entry to log file with given shot. No noise will be added, so noiseless
   * parameters will equal noisy ones.
   * Shot will not be executed, but noise information and agent names will be
   * written.
   */
  void write(const GameShot &gs, double duration = 0.0) {
    write(gs, gs.params, duration);
  }

  /**
   * Add Entry to log file with given shot. No noise will be added.
   * Shot will not be executed, but noise information and agent names will be
   * written.
   * \param gs Shot to write to file (actual shot executed).
   * \param noiselessParams Shot parameters given by AI before noise was added.
   * \param duration Simulation time in seconds for shot, or zero to avoid
   * including simulation time.
   */
  void write(const GameShot &gs, const ShotParams &noiselessParams,
             double duration = 0.0);

  /**
   * Add a brief state entry to log file. Entry will only include location of
   * balls of the table
   * and will not include any game-related information. Mostly useful for
   * physics experimentation and debugging. */
  void write(const TableState &ts);

  /**
   * Add a brief shot entry to log file. Entry will only include physical shot
   * parameters. Cue ball must already be positioned on the table.
   * Mostly useful for physics experimentation and debugging. */
  void write(const ShotParams &sp);

  /**
   * Adds a comment to the log file. Comment will be displayed in game viewer
   * with next shot.
   */
  void comment(const std::string text);

  /**
   * Execute specified shot and if physically possible, write previous state and
   * the actual shot to file.
   * Note that that given game state will be modified.
   * \param gst Game State before the shot, will be modified to reflect the
   * state after the shot.
   * \param gsh Shot to execute. Noise will not be added.
   * \param noiselessParams Shot parameters before noise was applied.
   */
  ShotResult LogAndExecute(GameState &gst, const GameShot &gsh,
                           const ShotParams &noiselessParams);

  /**
   * Add random noise, then execute specified shot and if physically possible,
   * write previous state and the actual shot to file.
   * Note that that given game state and shot will be modified.
   * \param gst Game State before the shot, will be modified to reflect the
   * state after the shot.
   * \param gsh Shot to execute. Noise will be added, and gsh will be modified
   * to reflect parameters with noise added.
   */
  ShotResult LogAndExecute(GameState &gst, GameShot &gsh);

  ~LogWriter() { delete _ofs; };

protected:
  /**
   * Utility function to write the GTYPE header to the stream
   */
  void writeHeader(const GameType gt);
  /**
   * Utility function to serialize a,b,theta,phi,v to the stream
   */
  void writeShotParams(const ShotParams &sp);

  std::string _agentName, _opponentName;
  const Noise *_noise;
  std::ofstream *_ofs;
  std::ostream &_os;
};
/*
   class LogReader {
// TODO: Specify and implement log reader
public:
//      LogReader(istream &is);
};*/
}

#endif
