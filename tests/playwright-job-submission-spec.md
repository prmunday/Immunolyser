# Immunolyser job submission via Playwright MCP — form spec

Notes from driving the `/initialiser` page end-to-end with the Playwright MCP
browser tools, for future debugging/regression-testing sessions. Update this
as the form changes or as more of the flow gets exercised.

Target used while writing this: test VM `http://118.138.243.174:5000`.

## Gotcha: full-page snapshot is too large

`browser_snapshot` on `/initialiser` returns 500k+ characters and fails with
a token-limit error. Use `browser_run_code_unsafe` with `page.$$eval(...)`
to query just the form elements instead, e.g.:

```js
async (page) => {
  const inputs = await page.$$eval('input, select, textarea', els => els.map(e => ({
    tag: e.tagName, type: e.type, id: e.id, name: e.name, checked: e.checked,
    value: (e.tagName==='SELECT'?undefined:e.value), visible: e.offsetParent !== null
  })).filter(e => e.visible));
  return JSON.stringify(inputs, null, 1);
}
```

## Step 1 — `/initialiser` page, sample/control setup

Fields present on initial page load:

| id | type | notes |
|---|---|---|
| `sample_name-1` | text | first sample's name |
| `sample_file-1` | file | first sample's replicate file(s) — accepts multiple files via `setInputFiles([...])` |
| `add_sample` (button) | — | click to append another `sample_name-N` / `sample_file-N` pair |
| `remove_sample` (button) | — | removes last sample row |
| `control_name-1` | text (hidden) | pre-filled `"Control"`, stays hidden — don't need to touch it |
| `control_file` | file | control/blank sample file(s) |
| `motif_length` | select | e.g. `8`, `9` |
| `species` | select | `Human` / `Mouse` |
| `mhc_class` | select | `I` / `II` (options depend on species) |
| `alleles_search` | text | filters the `alleles_list` multi-select live |
| `alleles_list` | select multiple | double-click an `<option>` to add it to the `alleles` textarea — **do not use `selectOption` first**, it leaves stray selections that get added alongside on the next dblclick |
| `alleles` | textarea | comma-separated allele names — **safest to just `page.fill('#alleles', 'H-2IAb')` directly** rather than relying on the list-double-click UI |
| `useFullDB` | select | "Use default list" / "Use full list" — only appears after species+class chosen in some flows; re-check visibility before relying on it |
| `submit` | input[type=submit] | value = `"Submit Job"` |

### IMPORTANT — one sample per replicate, not grouped

Immunolyser treats each uploaded "sample" as its own top-level output
directory. If the original job had 4 physically separate sample folders
(e.g. `Control_1_II`, `Control_2_II`, `Infected_1_II`, `Infected_2_II`), you
must recreate that as **4 separate `sample_name-N`/`sample_file-N` pairs**,
each with exactly one file — do NOT bundle multiple replicate files under one
sample name unless you've confirmed via server-side directory listing
(`ls /pvol/<taskId>/`) that's actually how the original samples were
structured. Check the directory names, not assumptions.

### Verified working sequence (mouse Class II example)

```js
await page.fill('#sample_name-1', 'Control_1_II');
await page.$('#sample_file-1').then(el => el.setInputFiles(['<path>/Control_1_II.csv']));
await page.fill('#sample_name-2', 'Control_2_II');
await page.$('#sample_file-2').then(el => el.setInputFiles(['<path>/Control_2_II.csv']));
await page.click('#add_sample'); // repeat per extra sample needed
await page.fill('#sample_name-3', 'Infected_1_II');
// ...
await page.$('#control_file').then(el => el.setInputFiles(['<path>/blank_control.csv']));

await page.selectOption('#motif_length', '9');
await page.selectOption('#species', 'Mouse');
await page.selectOption('#mhc_class', 'II');
await page.fill('#alleles_search', 'IAb');        // filter list to sanity-check the allele exists
await page.fill('#alleles', 'H-2IAb');            // set directly, don't rely on dblclick-to-add
```

## Step 2 — `/job-confirmation/<taskId>` page

Clicking `#submit` on the initialiser immediately creates the job (files are
uploaded/processed server-side) and redirects to
`http://<host>/job-confirmation/<taskId>` — the taskId is already assigned
at this point, before job name/email are even entered. Fields here:

| id | type | notes |
|---|---|---|
| `jobNameInput` | text | job's display name |
| `emailInput` | email | notification email (job completion / failure) |

```js
await page.fill('#jobNameInput', 'My Job Name');
await page.fill('#emailInput', 'someone@monash.edu');
```

There are **no prediction-tool checkboxes anywhere in the flow** — tools are
auto-selected server-side based on species + MHC class + chosen allele
compatibility (confirmed against `Immunolyser2.0_Allele_Dictionary.csv`).
E.g. Mouse + Class II + `H-2IAb` always runs `MixMHC2pred` + `NetMHCpanII`
since those are the only two tools marked compatible with mouse alleles in
the dictionary — matches what was observed in the original failed job's
celery logs (`Prediction tools selected: [MixMHC2pred, NetMHCpanII]`).

### Final step — attach job name/email to the already-running job

**Important:** the job starts processing immediately when `#submit` is
clicked on the initialiser page — `/job-confirmation/<taskId>` is not a
"review before running" step, it's a status/notification page for a job
that's already executing in the background. Job name/email are optional
metadata attached afterward via a `Submit Email` button (no stable id seen —
select by text: `button:has-text("Submit Email")`).

```js
await page.click('button:has-text("Submit Email")');
```

On success the page replaces its content with:
```
Request for Immunolyser report has been received. Task ID is <taskId>
...
Click here to refresh status.
Details submitted successfully!
```

There's also a `#checkStatusBtn` (text: "Click here to refresh status.") to
poll job status without a full page reload.

## Open questions / to fill in as discovered
- What `#checkStatusBtn` actually shows/does on click (in-place status
  update vs. navigation) — not yet driven.
- What the completed-job report page looks like / how to navigate to it
  once status is SUCCESS.
