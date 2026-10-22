# Submission order for Josh, from 16:00 EDT

All files named below are in `docs/submission/` unless a full path is given.

## 1. Go public

1. Confirm that https://github.com/joshbarron-36/altermagagent is public: open it in a private browser window while signed out.

## 2. LinkedIn post

1. Start a new post on LinkedIn.
2. Add the video: `altermagagent_final.mp4` if the voiced cut exists, else the silent cut `~/development/osi-hackathon/video/altermagagent_silent_v5.mp4` (the same file is in `docs/video/`).
3. Set visibility to "Anyone".
4. Paste the full text of `linkedin_post.txt`.
5. Post, then open the post and copy its URL (menu "..." then "Copy link to post").

## 3. Moment form, field by field

1. **Social media post URL:** paste the LinkedIn post URL from step 2.5.
2. **Write-up text box:** paste the whole of `writeup_moment.md`.
   Then replace each marker with the image, by dragging the file onto the marker line and deleting the marker text:
   - `[IMAGE 1]` with `images/01_cover_branch_graph.png`
   - `[IMAGE 2]` with `images/02_wannier_vs_vasp.png`
   - `[IMAGE 3]` with `images/03_mnte_bands.png`
   Keep the italic caption line under each image.
3. **Repository URL:** `https://github.com/joshbarron-36/altermagagent`, and tick "has a README".
4. **Cover and gallery images,** in this order (the first is the cover):
   1. `images/01_cover_branch_graph.png`
   2. `images/02_wannier_vs_vasp.png`
   3. `images/03_mnte_bands.png`
   4. `images/04_structure_moments.png`
   5. `images/05_marimo.png` (SKIP: this capture is from an older interface; do not upload)
5. **Tracks:**
   - Valency: yes. The literature module (`altermagagent/literature.py`) calls the Valency MCP endpoint https://labs.valency.io/mcp, and the chat has the `search_literature` tool.
   - Marimo: yes. The walkthrough notebook is `notebooks/altermagagent_walkthrough.py`.
   - Adaption: no. It is not used.
6. **Rules confirmation:** tick it.
7. **Review and submit:** read the preview, then submit.

## 4. After submitting

The entry stays editable until the deadline at 20:00 EDT, so a later edit can swap the video link (for example to a voiced cut) without a new entry.
