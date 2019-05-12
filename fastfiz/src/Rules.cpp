#include "Rules.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace Pool;

/*************************************************************************
 *************************************************************************
 **********BASE CLASS IMPLEMENTATION *************************************
 ************************************************************************/

/////////////////////////////////////////////
// Factory function
/////////////////////////////////////////////
GameState *GameState::Factory(istream &infoStream) {
  int gameType;
  infoStream >> gameType;
  GameType game_type((GameType)gameType); // Should cast things
  switch ((int)game_type) {
  case GT_EIGHTBALL: {
    EightBallState *newGame = new EightBallState(infoStream);
    return newGame;
    break;
  }

  case GT_NINEBALL: {
    NineBallState *newGame = new NineBallState(infoStream);
    return newGame;
    break;
  }

  default:
    cerr << "Unidentified gameType." << endl;
    cerr << "GameType: " << game_type << endl;
    throw "Unidentified gameType";
  }
  return NULL;
}

GameState *GameState::RackedState(GameType gameType) {
  switch ((int)gameType) {
  case GT_EIGHTBALL: {
    EightBallState *newGame = new EightBallState();
    return newGame;
  }

  default:
    cerr << "Unidentified gameType." << endl;
    break;
  }
  return NULL;
}

///////////////////////////////////////////////
// IO and Strings
///////////////////////////////////////////////
istream &Pool::operator>>(istream &aStream, GameState &gs) {
  gs.importFromStream(aStream);
  return aStream;
}

ostream &Pool::operator<<(ostream &stream, const GameState &gs) {
  gs.toStream(stream);
  return stream;
}

void GameState::importFromStream(istream &sourceStream) {
  sourceStream >> _turnType >> _timeLeft >> _timeLeftOpp >> _curPlayerStarted;
  _tableState.fromStream(sourceStream);
  if (!isLegalTurnType())
    throw "Illegal turn type in input!\n";
}

string GameState::toString() {
  ostringstream output;
  toStream(output);
  return output.str();
}

void GameState::toStream(ostream &out) const {
  out << gameType() << " " << _turnType << " " << _timeLeft << " ";
  out << _timeLeftOpp << " " << _curPlayerStarted << " ";
  _tableState.toStream(out);
}

//////////////////////////////////////////////////////////
// Shot Execution
//////////////////////////////////////////////////////////
ShotResult GameState::executeShot(const GameShot &shot, Shot **shotObj) {
  _switchedSides = false;
  cerr << "Trying to execute shot." << endl;

  if (shot.timeSpent && _timeLeft) {
    _timeLeft -= shot.timeSpent;
    if (_timeLeft < 0) {
      return SR_TIMEOUT;
    }
  }

  PreProcessCode code = preProcess(shot);
  cerr << "Past preProcess, value was: " << code << endl;

  if (code == PPC_BADPARAMS)
    return SR_BAD_PARAMS;
  if (code == PPC_G_NOEXECUTE)
    return _switchedSides ? SR_OK_LOST_TURN : SR_OK;
  if (code == PPC_G_PLACECUE) {
    cerr << "Placing ball at " << shot.cue_x << "," << shot.cue_y << endl;
    _tableState.setBall(Ball::CUE, Ball::STATIONARY, shot.cue_x, shot.cue_y);
  }
  // If there is a bad placement, that will be caught by tablestate's excute
  // shot
  // which will throw an exception.

  Shot *shotObject;
  try {
    shotObject = _tableState.executeShot(shot.params);
  } catch (BadShotException badshot) {
    // enum Type { OK, MISSED_BALL, INVALID_CUE_PARAMS, POOLFIZ_ERROR,
    // UNKNOWN_ERROR };
    cerr << "Got here." << endl;
    return SR_SHOT_IMPOSSIBLE;
  }

  if (shotObj) {
    *shotObj = shotObject;
  }

  const vector<Event *> &eventVec(shotObject->getEventList());
  try {
    processShot(eventVec, shot);
  } catch (BadShotException badshot) {
    if (!shotObj)
      delete shotObject;
    return SR_SHOT_IMPOSSIBLE;
  }

  if (!shotObj)
    delete shotObject;
  return _switchedSides ? SR_OK_LOST_TURN : SR_OK;
}

void GameState::switchSides() {
  int temp = _timeLeft;
  _timeLeft = _timeLeftOpp;
  _timeLeftOpp = temp;
  _curPlayerStarted = !(_curPlayerStarted);
  _switchedSides = !_switchedSides;
}

/*************************************************************************
 *************************************************************************
 **********EIGHT BALL IMPLEMENTATION *************************************
 ************************************************************************/

void EightBallState::importFromStream(istream &sourceStream) {
  // Call base class import
  GameState::importFromStream(sourceStream);

  // Import this class's parameters from end of string
  sourceStream >> _openTable >> _solids;
}

istream &Pool::operator>>(istream &aStream, EightBallState &gs) {
  gs.importFromStream(aStream);
  return aStream;
}

void EightBallState::switchSides() {
  GameState::switchSides();
  _solids = !_solids;
}

GameState::PreProcessCode EightBallState::preProcess(const GameShot &shot) {
  if (shot.decision == DEC_CONCEDE) {
    switchSides();
    _turnType = TT_WIN;
    return PPC_G_NOEXECUTE;
  }

  PreProcessCode returnValue = PPC_G_NORMAL;

  switch ((int)_turnType) {
  // Decision turn types
  case TT_EIGHTBALL_FOUL_ON_BREAK:
    switch ((int)shot.decision) {
    case DEC_EIGHTBALL_RERACK_OPP_SHOOT:
      switchSides();
    // Falling thorough...
    case DEC_RERACK:
      rack();
      _turnType = TT_BREAK;
      break;
    case DEC_KEEP_SHOOTING:
      _turnType = TT_BEHIND_LINE;
      break;
    default:
      return PPC_BADPARAMS;
    } // end shot.decision switch
    return PPC_G_NOEXECUTE;
  // No break needed, return takes care of it
  case TT_EIGHTBALL_8BALL_POCKETED_ON_BREAK:
    switch ((int)shot.decision) {
    case DEC_RERACK:
      rack();
      _turnType = TT_BREAK;
      break;
    case DEC_KEEP_SHOOTING:
      _turnType = (_tableState.getBall(Ball::CUE).isInPlay() ? TT_NORMAL
                                                             : TT_BEHIND_LINE);
      break;
    default:
      return PPC_BADPARAMS;
    } // end shot.decision switch
    return PPC_G_NOEXECUTE;

  // Shot turn types
  case TT_BREAK:
  case TT_BEHIND_LINE: {
    double headstring = _tableState.getTable().getHeadString();
    if (shot.cue_y < headstring)
      return PPC_BADPARAMS;
  }
  // Intentionally falling through...
  case TT_BALL_IN_HAND:
    returnValue = GameState::PPC_G_PLACECUE;
  // Intentionally falling through
  case TT_NORMAL:
    if (_turnType != TT_BREAK) {
      if (shot.ball > Ball::UNKNOWN_ID || shot.ball < Ball::CUE)
        return PPC_BADPARAMS;
      if (shot.pocket > Table::UNKNOWN_POCKET || shot.pocket < Table::SW)
        return PPC_BADPARAMS;
    }
    return returnValue;

  default:
    throw "No case in preProcess of EightBallState";
    break;

  } // end switch _turnType

  return returnValue;
}

void EightBallState::processShot(const vector<Event *> &eventList,
                                 const GameShot &gameShot) {
  vector<Event *>::const_iterator itr;
  itr = eventList.begin();

  // Special case out the first event: should be CUE_STRIKE, and not a MISCUE
  if ((*itr)->getType() != Event::CUE_STRIKE) {
    throw BadShotException(BadShotException::MISSED_BALL);
  }
  itr++;

  switch ((int)_turnType) {
  case TT_BREAK: {
    cerr << "In turn-type BREAK!" << endl;
    int rail_collisions = 0;
    int pocketed_balls = 0;
    bool eight_ball_pocketed = false;
    bool scratch = false;
    for (; itr != eventList.end(); itr++) {
      if ((*itr)->getType() == Event::POCKETED) {
        if ((*itr)->getBall1() == Ball::EIGHT)
          eight_ball_pocketed = true;
        if ((*itr)->getBall1() == Ball::CUE)
          scratch = true;
        else
          pocketed_balls++;
      } else if ((*itr)->getType() == Event::RAIL_COLLISION)
        rail_collisions++; // Should check for different balls & note cue ball
    }                      // end for loop

    bool legal_shot = false;

    if (!scratch && !eight_ball_pocketed &&
        (rail_collisions > 3 || pocketed_balls > 0))
      legal_shot = true;
    if (legal_shot) {
      _turnType = TT_NORMAL;
      _openTable = true;
    } else if (!scratch && !eight_ball_pocketed && !legal_shot) {
      _turnType = TT_EIGHTBALL_FOUL_ON_BREAK; // Options for decision: Normal /
                                              // Re-Rack / re-rack opp break
      _openTable = true;
      switchSides();
    } else if (scratch && !eight_ball_pocketed) {
      _openTable = true;
      _turnType = TT_BEHIND_LINE;
      switchSides();
    } else if (eight_ball_pocketed) {
      if (scratch)
        switchSides();
      _turnType = TT_EIGHTBALL_8BALL_POCKETED_ON_BREAK;
      _openTable = true;
      _tableState.spotBall(Ball::EIGHT);
    };
    // DEBUG CODE - never lose turn after break
    /*if (_turnType != TT_NORMAL) {
      switchSides();
      rack();
      _turnType = TT_BREAK;
      }*/
  } // End of case TT_BREAK
  break;

  case TT_NORMAL:
  case TT_BALL_IN_HAND:
  case TT_BEHIND_LINE:
    cerr << "WE've gotten to the TT_NORMAL, BallInHand,BehindLine cases"
         << endl;
    {
      bool pocketed_legal_ball =
          false;                   // Set to true if a legal ball was pocketed.
      bool rail_collision = false; // Set to true if a ball hits a rail after
                                   // the cue contacs the first object ball
      bool scratch = false;        // set to true if the cue ball is sunk.
      bool eight_ball_pocketed =
          false; // set to true if the eight ball is pocketed
      bool first_ball_collision =
          false; // set to true after the first ball collision occurs
      // TODO: Check use of safety_shot
      // bool safety_shot = false;
      // if (gameShot.ball == Ball::UNKNOWN_ID) safety_shot=true;
      bool foul = false; // For any other sorts of fouls
      bool past_head_string = false;
      Ball::Type first_ball = Ball::CUE;
      bool all_sunk = true;
      if (!_openTable) {
        int name = _solids ? Ball::ONE : Ball::NINE;
        int end_name = _solids ? Ball::EIGHT : Ball::UNKNOWN_ID;
        for (; name < end_name; name++) {
          if (!_tableState.getBall(Ball::Type(name)).isPocketed())
            all_sunk = false;
        }
      } else {
        all_sunk = false;
      }

      for (; itr != eventList.end(); itr++) {
        if (_turnType == TT_BEHIND_LINE && !first_ball_collision &&
            !past_head_string) {
          Ball ball1 = (*itr)->getBall1Data();
          Point ball1Location = ball1.getPos();
          double y_val = ball1Location.y;
          double headstring = _tableState.getTable().getHeadString();
          if (y_val > headstring)
            past_head_string = true;
        }
        switch ((int)(*itr)->getType()) {
        case Event::POCKETED: {
          Ball::Type ball1 = (*itr)->getBall1();
          Table::Pocket pocket;
          PocketedEvent *pock_event;
          pock_event = dynamic_cast<PocketedEvent *>(*itr);
          pocket = pock_event->getPocket();
          switch ((int)ball1) {
          case Ball::CUE:
            scratch = true;
            break;

          case Ball::ONE:
          case Ball::TWO:
          case Ball::THREE:
          case Ball::FOUR:
          case Ball::FIVE:
          case Ball::SIX:
          case Ball::SEVEN:
            if (_solids)
              all_sunk = false;
            if (gameShot.ball == ball1 && pocket == gameShot.pocket &&
                (_openTable || _solids)) {
              pocketed_legal_ball = true;
            }
            break;

          case Ball::EIGHT:
            eight_ball_pocketed = true;
            if (all_sunk && gameShot.ball == ball1 &&
                pocket == gameShot.pocket) {
              pocketed_legal_ball = true;
            }
            break;

          case Ball::NINE:
          case Ball::TEN:
          case Ball::ELEVEN:
          case Ball::TWELVE:
          case Ball::THIRTEEN:
          case Ball::FOURTEEN:
          case Ball::FIFTEEN:
            if (!_solids)
              all_sunk = false;
            if (gameShot.ball == ball1 && pocket == gameShot.pocket &&
                (_openTable || !_solids)) {
              pocketed_legal_ball = true;
            }
            break;

          } // switch on ball1
        }   // case Event::POCKETED
        break;

        case Event::RAIL_COLLISION:
          if (first_ball_collision)
            rail_collision = true;
          break;

        case Event::BALL_COLLISION:
          if (!first_ball_collision) {
            first_ball_collision = true;
            first_ball = (*itr)->getBall2();
          }
          break;

        } // end Event switch

      } // end events  for loop

      // Evaluate the boolean values to foul, scratch, and safety_shot
      if (!first_ball_collision)
        foul = true;
      if (!rail_collision && !pocketed_legal_ball)
        foul = true;
      if (_turnType == TT_BEHIND_LINE && !past_head_string)
        foul = true;
      if (!_openTable && _solids && !all_sunk &&
          !(first_ball > Ball::CUE && first_ball < Ball::EIGHT))
        foul = true;
      if (!_openTable && !_solids && !all_sunk &&
          !(first_ball > Ball::EIGHT && first_ball < Ball::UNKNOWN_ID))
        foul = true;
      if (all_sunk && first_ball != Ball::EIGHT)
        foul = true;
      cerr << "Open table: " << _openTable << "Solids: " << _solids
           << "All sunk " << all_sunk << "First ball (enum): " << first_ball
           << "EightBallPocketed " << eight_ball_pocketed << endl;

      if (eight_ball_pocketed) {
        _turnType = TT_WIN;
        if (!all_sunk || scratch || foul)
          switchSides();
      } else if (scratch || foul) {
        switchSides();
        _turnType = TT_BALL_IN_HAND;
        return;
      } else if (pocketed_legal_ball) {
        _turnType = TT_NORMAL;
        if (_openTable) {
          _openTable = false;
          _solids = gameShot.ball < Ball::EIGHT;
        }
      } else {
        _turnType = TT_NORMAL;
        switchSides();
      }

    } // end TT_NORMAL and TT_BALL_IN_HAND and TT_BEHIND_LINE
    break;

  default:
    throw "No _gameType in processShot";
    break;

  } // switch _gameType

} // EightBallState::processShot()

void EightBallState::rack() {
  const Table &t = _tableState.getTable();
  double w = t.getWidth();
  double x = t.getFootSpot().x;
  double y = t.getFootSpot().y;
  double r = Ball::BALL_RADIUS;
  double dither = 0.00005;

  double dy = r * sqrt(3) * (1 + dither);
  double dx = r * (1 + dither);
  _tableState.setBall(Ball::CUE, Ball::STATIONARY, w / 2, t.getHeadString());
  _tableState.setBall(Ball::ONE, Ball::STATIONARY, x, y);
  //_tableState.getBall(Ball::ONE).setState(Ball::SLIDING);
  //_tableState.getBall(Ball::ONE).setSpin(Vector(180,0,0));

  _tableState.setBall(Ball::TWO, Ball::STATIONARY, x - dx, y - dy);
  _tableState.setBall(Ball::THREE, Ball::STATIONARY, x + dx, y - dy);
  _tableState.setBall(Ball::FOUR, Ball::STATIONARY, x - dx * 2, y - dy * 2);
  _tableState.setBall(Ball::FIVE, Ball::STATIONARY, x, y - dy * 2);
  _tableState.setBall(Ball::SIX, Ball::STATIONARY, x + dx * 2, y - dy * 2);
  _tableState.setBall(Ball::SEVEN, Ball::STATIONARY, x - dx * 3, y - dy * 3);
  _tableState.setBall(Ball::EIGHT, Ball::STATIONARY, x - dx, y - dy * 3);
  _tableState.setBall(Ball::NINE, Ball::STATIONARY, x + dx, y - dy * 3);
  _tableState.setBall(Ball::TEN, Ball::STATIONARY, x + dx * 3, y - dy * 3);
  _tableState.setBall(Ball::ELEVEN, Ball::STATIONARY, x - dx * 4, y - dy * 4);
  _tableState.setBall(Ball::TWELVE, Ball::STATIONARY, x - dx * 2, y - dy * 4);
  _tableState.setBall(Ball::THIRTEEN, Ball::STATIONARY, x, y - dy * 4);
  _tableState.setBall(Ball::FOURTEEN, Ball::STATIONARY, x + dx * 2, y - dy * 4);
  _tableState.setBall(Ball::FIFTEEN, Ball::STATIONARY, x + dx * 4, y - dy * 4);

  _tableState.addNoise(dither);
}

namespace Pool {
  string getRulesVersion() {
    return "Rules version 0.1 built " __DATE__ " " __TIME__;
  }
}

void EightBallState::toStream(ostream &out) const {
  GameState::toStream(out);
  out << " " << _openTable << " " << _solids;
}

GameState *EightBallState::clone() {
  EightBallState *gs = new EightBallState(*this);
  return gs;
}

/*************************************************************************
 *************************************************************************
 **********NINE BALL IMPLEMENTATION *************************************
 ************************************************************************/

void NineBallState::importFromStream(istream &sourceStream) {
  // Call base class import
  GameState::importFromStream(sourceStream);

  // Import this class's parameters from end of string
  sourceStream >> _fouls >> _opp_fouls;
}

istream &Pool::operator>>(istream &aStream, NineBallState &gs) {
  gs.importFromStream(aStream);
  return aStream;
}

void NineBallState::switchSides() {
  GameState::switchSides();
  int tmp = _opp_fouls;
  _opp_fouls = _fouls;
  _fouls = tmp;
}

GameState::PreProcessCode NineBallState::preProcess(const GameShot &shot) {
  if (shot.decision == DEC_CONCEDE) {
    switchSides();
    _turnType = TT_WIN;
    return PPC_G_NOEXECUTE;
  }

  PreProcessCode returnValue = PPC_G_NORMAL;

  switch ((int)_turnType) {
  // Decision turn types
  case TT_NINEBALL_PUSH_OUT:
    if (shot.decision == DEC_KEEP_SHOOTING) {
      _turnType = TT_NORMAL;
      return PPC_G_NOEXECUTE;
    } else if (shot.decision == DEC_NINEBALL_PUSH_OUT) {
      _turnType = TT_NORMAL;
      switchSides();
      return PPC_G_NOEXECUTE;
    } else {
      return PPC_BADPARAMS;
    }

  // Shot turn types
  case TT_BREAK: {
    double headstring = _tableState.getTable().getHeadString();
    if (shot.cue_y < headstring)
      return PPC_BADPARAMS;
  }
  // Intentionally falling through...
  case TT_BALL_IN_HAND:
    returnValue = GameState::PPC_G_PLACECUE;
  // Intentionally falling through
  case TT_NORMAL:
    if (shot.decision != DEC_NO_DECISION)
      return PPC_BADPARAMS;
  case TT_NINEBALL_FIRST_SHOT:
    return returnValue;

  default:
    throw "No case in preProcess of NineBallState";
    break;

  } // end switch _turnType

  return returnValue;
}

void NineBallState::processShot(const vector<Event *> &eventList,
                                const GameShot &gameShot) {
  vector<Event *>::const_iterator itr;
  itr = eventList.begin();

  // Special case out the first event: should be CUE_STRIKE, and not a MISCUE
  if ((*itr)->getType() != Event::CUE_STRIKE) {
    throw BadShotException(BadShotException::MISSED_BALL);
  }
  itr++;

  Ball::Type lowestBall = Ball::NINE;

  for (Ball::Type b = Ball::ONE; b < Ball::NINE; b = (Ball::Type)((int)b + 1)) {
    if (_tableState.getBall(b).isInPlay()) {
      lowestBall = b;
      break;
    }
  }

  int rail_collisions = 0;
  int pocketed_balls = 0;
  bool nine_ball_pocketed = false;
  bool scratch = false;
  Ball::Type first_ball_hit = Ball::UNKNOWN_ID;
  for (; itr != eventList.end(); itr++) {
    if ((*itr)->getType() == Event::POCKETED) {
      Ball::Type b = (*itr)->getBall1();
      if (b == Ball::NINE)
        nine_ball_pocketed = true;
      if (b == Ball::CUE)
        scratch = true;
      else {
        if (b < lowestBall)
          lowestBall = b;
        pocketed_balls++;
      }
    } else if ((*itr)->getType() == Event::RAIL_COLLISION) {
      Ball::Type b = (*itr)->getBall1();
      if ((_turnType == TT_BREAK && b != Ball::CUE) ||
          (_turnType != TT_BREAK && first_ball_hit != Ball::UNKNOWN_ID)) {
        rail_collisions++;
      }
    } else if ((*itr)->getType() == Event::BALL_COLLISION) {
      Ball::Type b = (*itr)->getBall2();
      if (first_ball_hit == Ball::UNKNOWN_ID) {
        first_ball_hit = b;
      }
    }
  } // end for loop
  bool foul = scratch;
  if (gameShot.decision == DEC_NINEBALL_PUSH_OUT) {
    pocketed_balls = 0; // pocketed balls don't count.
  } else {
    foul = foul || (first_ball_hit != lowestBall);
    foul = foul || ((rail_collisions < (_turnType == TT_BREAK ? 4 : 1)) &&
                    !pocketed_balls);
  }

  if (!foul)
    _fouls = 0; // Reset consecutive foul counter

  bool keepTurn = pocketed_balls && !foul;
  if (!keepTurn)
    switchSides();

  if (nine_ball_pocketed && !keepTurn) {
    _tableState.spotBall(Ball::NINE);
    nine_ball_pocketed = false;
  }

  if (foul) {
    if (++_opp_fouls == 3)
      _turnType = TT_WIN;
    else
      _turnType = TT_BALL_IN_HAND;
  } else if (nine_ball_pocketed) {
    _turnType = TT_WIN;
  } else if (_turnType == TT_BREAK) {
    _turnType = TT_NINEBALL_FIRST_SHOT;
  } else if (gameShot.decision == DEC_NINEBALL_PUSH_OUT) {
    _turnType = TT_NINEBALL_PUSH_OUT;
  } else {
    _turnType = TT_NORMAL;
  }
} // NineBallState::processShot()

void NineBallState::rack() {
  const Table &t = _tableState.getTable();
  double w = t.getWidth();
  double x = t.getFootSpot().x;
  double y = t.getFootSpot().y;
  double r = Ball::BALL_RADIUS;
  double dither = 0.00005;

  double dy = r * sqrt(3) * (1 + dither);
  double dx = r * (1 + dither);
  _tableState.setBall(Ball::CUE, Ball::STATIONARY, w / 2, t.getHeadString());
  _tableState.setBall(Ball::ONE, Ball::STATIONARY, x, y);
  _tableState.setBall(Ball::TWO, Ball::STATIONARY, x - dx, y - dy);
  _tableState.setBall(Ball::THREE, Ball::STATIONARY, x + dx, y - dy);
  _tableState.setBall(Ball::FOUR, Ball::STATIONARY, x - dx * 2, y - dy * 2);
  _tableState.setBall(Ball::NINE, Ball::STATIONARY, x, y - dy * 2);
  _tableState.setBall(Ball::FIVE, Ball::STATIONARY, x + dx * 2, y - dy * 2);
  _tableState.setBall(Ball::SIX, Ball::STATIONARY, x + dx, y - dy * 3);
  _tableState.setBall(Ball::SEVEN, Ball::STATIONARY, x - dx, y - dy * 3);
  _tableState.setBall(Ball::EIGHT, Ball::STATIONARY, x, y - dy * 4);

  _tableState.addNoise(dither);
}

void NineBallState::toStream(ostream &out) const {
  GameState::toStream(out);
  out << " " << _fouls << " " << _opp_fouls;
}

GameState *NineBallState::clone() {
  NineBallState *gs = new NineBallState(*this);
  return gs;
}
