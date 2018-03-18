#include "AIBase.h"

using namespace std;
using namespace Pool;

GameShot &Pool::AIBase::computeShot(const GameState &gs,
                                    const Pool::Noise *const p_noise) {
  _stopwatch->restart();
  gameState = &gs;
  noise = p_noise;
  shot.decision = DEC_NO_DECISION;
  shot.pocket = Table::UNKNOWN_POCKET;
  shot.ball = Ball::UNKNOWN_ID;
  retries = MAX_RETRIES;
  // cerr << "Computing shot for game state " << gs << endl;
  if (!forGame(gs.gameType())) {
    cerr << "Can't handle this game, conceding!" << endl;
    shot.decision = DEC_CONCEDE;
  } else if (gs.shotRequired()) {
    shoot();
  } else if (gs.decisionAllowed()) {
    shot.decision = decide();
  } else if (gs.isTerminal()) {
    cerr << "Shot requested for terminal state!" << endl;
    throw "Shot requested for terminal state!";
  }
  shot.timeSpent = _stopwatch->getElapsed();
  return shot;
}

GameShot &Pool::AIBase::reComputeShot() {
  if (retries--) {
    return shot;
  } else {
    shot.decision = DEC_CONCEDE;
    return shot;
  }
}

Decision Pool::AIBase::decide() {
  switch (gameState->getTurnType()) {
  case TT_EIGHTBALL_8BALL_POCKETED_ON_BREAK:
  case TT_EIGHTBALL_FOUL_ON_BREAK:
    return DEC_RERACK;
  default:
    return DEC_CONCEDE;
  }
}

void Pool::AIBase::shoot() {
  switch (gameState->getTurnType()) {
  case TT_BREAK:
    // cerr << "    Shooting break shot." << endl;
    breakShot();
    break;
  default:
    // cerr << "    Shooting non-break shot." << endl;
    otherShot();
    break;
  }
}
