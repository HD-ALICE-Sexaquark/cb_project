### Visualize only interesting events

In ``EndOfEventAction``:

```
auto eventManager = G4EventManager::GetEventManager();
  if (<some condition>) eventManager->KeepTheCurrentEvent();
```

In GUI:

```
/vis/disable
/run/beamOn 1000
/vis/enable
/vis/reviewKeptEvents
```
In the latter, repeatetly type ``cont``
