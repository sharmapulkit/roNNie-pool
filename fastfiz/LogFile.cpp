#include <iostream>
#include <sstream>
#include <string>
#include "LogFile.h"

using namespace std;
using namespace Pool;

/* TODO: Implement log writer */


void LogWriter::write(const GameState & gs)
{
    _os << "STATE " << gs << endl;
}

void LogWriter::writeShotParams(const ShotParams & sp) {
    _os << sp.a << " " << sp.b << " " << sp.theta << " " << sp.phi << " " << sp.v;
}

void LogWriter::write(const GameShot & gs, const ShotParams & noiselessParams, double duration)
{
    NoNoise nn;
    _os << "SHOT ";
    writeShotParams(gs.params);
    _os << " ";
    writeShotParams(noiselessParams);
    _os << " " << gs.ball << " " << gs.pocket << " " << gs.decision
        << " " << gs.cue_x << " " << gs.cue_y
        << " " << gs.timeSpent << " " << duration << " " ;
    if (_noise) 
        _noise->toStream(_os);
    else
        nn.toStream(_os);
    _os << " \"" << _agentName << "\" \"" << _opponentName << "\"" << endl;
}

void LogWriter::write(const TableState & ts)
{
    _os << "TSTATE ";
    ts.toStream(_os);
    _os << endl;
}

void LogWriter::write(const ShotParams & sp)
{
    _os << "TSHOT ";
    writeShotParams(sp);
    _os << endl;
}

void LogWriter::comment(const string text)
{
    _os << "COMMENT " << text << endl;
}

ShotResult LogWriter::LogAndExecute(GameState & gst, const GameShot & gsh, const ShotParams & noiselessParams)
{
    GameState *prevState = gst.clone();
    Shot* shotObj;
    ShotResult res = gst.executeShot(gsh,&shotObj);
    if (res == SR_OK || res == SR_OK_LOST_TURN) {
        write(*prevState);
        write(gsh,noiselessParams,shotObj->getDuration());
    }
    delete shotObj;
    delete prevState;
    return res;
}

ShotResult LogWriter::LogAndExecute(GameState & gst, GameShot & gsh)
{
    ShotParams nlp=gsh.params;
    if (_noise) {
        _noise->applyNoise(gsh.params);
    }
    return LogAndExecute(gst,gsh,nlp);
}

void LogWriter::writeHeader(const GameType gt)
{
    _os << "GTYPE " << gt << endl;
}

