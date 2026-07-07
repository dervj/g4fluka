#include "EventAction.hh"
#include "G4Event.hh"

EventAction::EventAction() {}
EventAction::~EventAction() {}

void EventAction::BeginOfEventAction(const G4Event*) {}

void EventAction::EndOfEventAction(const G4Event* event)
{
    // Progress printout every 50k events
    G4int evtNo = event->GetEventID();
    if (evtNo > 0 && evtNo % 50000 == 0)
        G4cout << ">>> Event " << evtNo << G4endl;
}
