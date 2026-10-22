# How the silent screen recording was produced

The silent cut of the video follows `docs/video_script.md` (version 3) scene by scene.
It was made by script on Josh's laptop on 22 October 2026; the video file is not in the repository.

## Output

- Video: `~/development/osi-hackathon/video/altermagagent_silent_v1.mp4`, 1920 by 1080, H.264, yuv420p, 25 frames per second, 3 minutes 58 seconds.
- Scene list: `~/development/osi-hackathon/video/scenes.json` gives each scene's name, title, start, duration and narration paragraph, for aligning Diana's voice.
- Each scene carries a caption band at the bottom with its title, so the silent cut reads on its own.
- The optional Wannier comparison scene is not in this cut, because Diana's comparison script has not landed.

## Steps

1. Start the interface with the hosted chat and the literature search: `source ~/.altermagagent/env`, then `.venv/bin/altermag gui` on port 8080.
2. Start headless Chromium at 1920 by 1080 with `--remote-debugging-port=9222 --hide-scrollbars` on http://localhost:8080/?project=CrSe_real.
3. Run `video/tools/record.mjs` with Node 24: it injects a drawn cursor, scrolls the storyline panel and clicks by text or by `data-branch` on the graph nodes.
4. For each browser scene it records DevTools screencast frames, sampled about five times per second, into `video/frames/<scene>/`.
5. Run `video/tools/terminal.py` to draw the terminal scenes (folder tree, import, accuracy, title card) as 1920 by 1080 frames with PIL.
6. Run `video/tools/compose.py`: it holds the newest frame at each 1/25 second, adds the caption band, encodes one clip per scene and concatenates them.

## What the frames show

- The index was not re-imported for the recording.
- The import scene shows the import command and the scan lines stored in the index from the last import: CrSe_real 13 runs, 1127 files, 3.46 s; MnTe 4 runs, 175 files, 1.3 s.
- The accuracy scene shows the real output of `python tests/accuracy.py`: phase 137 / 137, status 137 / 137, parent 135 / 137.
- The band scene cuts the wait while the server draws the band figure from vasprun.xml (about 40 to 90 seconds on this laptop).
- The 3D structure view opens in a new tab and needs internet, so the cut shows the static structure figure in the branch header instead of a rotation.
- The chat answer comes live from the hosted model; its length and wording differ from run to run.

## Rebuilding

The tools live outside the repository in `~/development/osi-hackathon/video/tools/`.
Run `node record.mjs`, then `terminal.py` and `compose.py` with the repository's `.venv/bin/python -I`; one scene can be re-recorded with `node record.mjs band,chat`.

## Version 3

Version 3 (3 minutes 26 seconds) re-records only summary, suggestions, msg (now Symmetries and convergence) and the new wannier_match scene with `video/tools/record_v3.mjs`.
`video/tools/compose_v3.py` encodes those four scenes and reuses the version 2 clips for the others; the cut and `scenes_v3.json` are in `docs/video/`.

Version 4 (3 minutes 26 seconds) re-records every interface scene with home paths shown as "~" or relative to the project, and redraws the terminal scenes with `terminal_v4.py`; `compose_v4.py` encodes it.
