# Project data handling

Read and follow [DATA_USE_RULES.md](DATA_USE_RULES.md) before working with data or preparing deliverables. These rules also apply to subagents, tool arguments and outputs, debug prints, figures, error messages, ZIP archives and git operations.

Keep clinical records and single-patient derived information in `restricted/` or existing `_restricted/` directories. Never return their contents to an AI conversation. Patient counts, rather than pair counts alone, determine whether a clinical summary can be shown.

Before any package, upload, commit or push, run `tools/small_cell_gate.py` on the actual payload and record its result. A failure blocks that transfer. Do not bypass hooks or exemptions by renaming fields or files. Use the guarded packager for new archives.

Use Chinese for author-facing reports and English for manuscript text. Preserve scientific estimates when changing disclosure layers. End each report with the R4 data-use self-check and the delivery SHA manifest.
