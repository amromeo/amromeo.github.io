// Render only recorded fits. Statistical fitting is in R; see tools/build-demo.R.
const byId = (id) => document.getElementById(id);
const form = byId("design");
const plot = byId("plot");
const tru = [24.120108046, 35.879891954];
const settings = {
  1: "Analyte-altering disease prevalence is set to 5% in clinic and 10% on the ward. Disease shifts results upward.",
  2: "Analyte-altering disease prevalence is set to 20% in clinic and 75% on the ward. Disease-related results overlap the upper end of the reference distribution.",
  3: "Analyte-altering disease prevalence is set to 20% in clinic and 75% on the ward. Ward disease can shorten or prolong the result.",
  4: "Analyte-altering disease prevalence is set to 40% in both settings. Setting provides no systematic difference in the proportion of reference-eligible patients."
};
const methods = {
  refineR_all: ["refineR · all results", "#586e9b"],
  refineR_clinic: ["refineR · outpatients", "#995022"],
  anchored: ["anchoredRI", "#006b60"],
  anchor_only: ["Screened alone", "#697774"],
  screened_120: ["Screened alone · 120", "#8c682f"],
  screened_400: ["Screened alone · 400", "#785280"],
  joint: ["Settings alone", "#41605e"],
  pooled_mixture: ["Pooled mixture", "#655777"],
  anchored_two: ["With screened patients · two components", "#367a75"]
};
const referenceMethods = ["screened_120", "screened_400"];
let rows = [];
let display = "limits";

function node(tag, cls, value) {
  const el = document.createElement(tag);
  if (cls) el.className = cls;
  if (value !== undefined) el.textContent = value;
  return el;
}

function parseFits(csv) {
  const lines = csv.trim().split(/\r?\n/);
  const fields = lines.shift().split(",");
  const data = lines.map(line => Object.fromEntries(line.split(",").map((v, i) => {
    const key = fields[i];
    if (key === "method" || key === "selection") return [key, v];
    if (key === "at_bound") return [key, v === "TRUE"];
    return [key, v === "" ? null : Number(v)];
  })));
  if (data.length !== 3648 || data.some(r => !methods[r.method] || !Number.isFinite(r.lower)
      || !Number.isFinite(r.upper) || r.lower >= r.upper)) {
    throw new Error("The fitted-results file is incomplete.");
  }
  return data;
}

function average(values) {
  return values.reduce((sum, v) => sum + v, 0) / values.length;
}

function summarize(sample) {
  if (!sample.length) return null;
  return {
    lower: average(sample.map(r => r.lower)),
    upper: average(sample.map(r => r.upper)),
    worst: average(sample.map(r => Math.max(Math.abs(r.lower - tru[0]), Math.abs(r.upper - tru[1])))),
    atBound: sample.filter(r => r.at_bound).length,
    count: sample.length
  };
}

function axisRange(intervals, error) {
  const values = intervals.filter(Boolean);
  if (error) {
    const largest = Math.max(1, ...values.map(r => r.worst));
    const step = largest > 12 ? 5 : largest > 5 ? 2 : largest > 2.5 ? 1 : .5;
    return [0, step * Math.max(3, Math.ceil(largest / step)), step];
  }
  const low = Math.min(22, ...values.map(r => r.lower));
  const high = Math.max(38, ...values.map(r => r.upper));
  const step = high - low > 30 ? 10 : 5;
  return [Math.floor(low / step) * step, Math.ceil(high / step) * step, step];
}

function makeRow(label, color, r, scale, reference = false) {
  const [min, max] = scale;
  const percent = (x) => (x - min) / (max - min) * 100;
  const row = node("div", reference ? "row reference-row" : "row");
  row.style.setProperty("--series", color);
  const heading = node("div", "method-name", label);
  const track = node("div", "track");
  track.setAttribute("aria-hidden", "true");
  const metric = node("div", "metric");
  if (!r) {
    heading.append(node("span", "row-value", "Needs screened results"));
    metric.textContent = "—";
  } else {
    heading.append(node("span", "row-value", r.lower.toFixed(2) + " – " + r.upper.toFixed(2) + " s"));
    if (display === "limits") {
      const band = node("span", "band");
      band.style.left = percent(tru[0]) + "%";
      band.style.width = percent(tru[1]) - percent(tru[0]) + "%";
      const estimate = node("span", "estimate");
      estimate.style.left = percent(r.lower) + "%";
      estimate.style.width = percent(r.upper) - percent(r.lower) + "%";
      track.append(band, estimate);
    } else if (!reference) {
      const bar = node("span", "errorbar");
      bar.style.width = percent(r.worst) + "%";
      track.append(bar);
    }
    if (reference) {
      metric.textContent = "known";
    } else {
      metric.append(document.createTextNode(r.worst.toFixed(2) + " "), node("small", "", "s"));
    }
  }
  row.append(heading, track, metric);
  return row;
}

function render() {
  const fields = new FormData(form);
  const scenario = Number(byId("scenario").value);
  const n = Number(fields.get("routine"));
  const anchors = Number(fields.get("anchors"));
  const rule = byId("selection").value;
  const draw = Number(byId("replicate").value);
  byId("scenario-note").textContent = settings[scenario];
  byId("sizes").textContent = (n / 4).toLocaleString("en-US") + " outpatients · " + (n * .75).toLocaleString("en-US") + " ward patients";
  byId("selection-note").textContent = rule === "clinical"
    ? "Reviewers are blinded to the target result for all screened cohorts."
    : "All screened cohorts retain only values from 26 through 34 seconds.";
  byId("result-title").textContent = display === "limits" ? "Estimated reference limits" : "Error in the reference limits";
  byId("plot-note").textContent = display === "limits"
    ? (draw ? "One simulated dataset. " : "Average limits across 12 datasets. ") + "The shaded band marks the known reference interval."
    : (draw ? "The bar shows the larger error from the lower and upper limits. " : "Each bar averages the larger limit error from each of 12 datasets. ") + "Smaller values mean more accurate limits.";
  byId("error-definition").textContent = draw
    ? "Compare each estimated limit with its known value. The error shown is the larger of the two absolute differences."
    : "For each of the 12 datasets, take the larger error from the lower and upper limits. The error column shows the average of those 12 numbers.";
  if (!rows.length) return;
  const selected = rows.filter(r => r.scenario === scenario && r.n_routine === n && (!draw || draw === r.replicate)
    && (r.selection === "any" || (r.selection === rule
      && (r.n_anchor === anchors || referenceMethods.includes(r.method)))));
  const active = ["refineR_all", "refineR_clinic", "anchored"];
  if (anchors > 0 && anchors !== 120 && anchors !== 400) active.push("anchor_only");
  active.push(...referenceMethods);
  if (byId("extra").checked) active.push("joint", "pooled_mixture", "anchored_two");
  const stats = active.map(id => summarize(selected.filter(r => r.method === id)));
  const scale = axisRange(stats, display === "error");
  const [min, max, step] = scale;
  const titleRow = node("div", "axis-row");
  const axis = node("div", "axis");
  axis.setAttribute("aria-hidden", "true");
  for (let v = min; v <= max + 1e-9; v += step) {
    const tick = node("span", "tick", String(v));
    tick.style.left = ((v - min) / (max - min) * 100) + "%";
    axis.append(tick);
  }
  const errorHeading = node("a", "axis-last", draw ? "Larger limit error" : "Average larger limit error");
  errorHeading.href = "help.html#limit-error";
  errorHeading.title = "How the limit error is calculated";
  titleRow.append(node("span", "", "Method / limits (s)"), axis, errorHeading);
  const items = [titleRow];
  if (display === "limits") items.push(makeRow("Known reference interval", "#32795e", {lower:tru[0], upper:tru[1]}, scale, true));
  active.forEach((id, i) => {
    const [label, color] = methods[id];
    const name = id === "anchor_only" ? label + " · " + anchors : label;
    items.push(makeRow(name, color, stats[i], scale));
  });
  plot.replaceChildren(...items);
  const warnings = active.flatMap((id, i) => stats[i]?.atBound ? [methods[id][0] + " (" + stats[i].atBound + "/" + stats[i].count + ")"] : []);
  byId("diagnostics").textContent = warnings.length ? "A numerical fitting limit was reached for " + warnings.join(", ") + ". See Help for the meaning of this warning." : "";
  byId("announcement").textContent = "Updated " + (draw ? "dataset " + draw : "12-dataset average") + " for " + n + " routine and " + anchors + " screened patients.";
}

async function loadFits() {
  try {
    const response = await fetch("assets/fits.csv?v=screened120400");
    if (!response.ok) throw new Error("Results are unavailable.");
    rows = parseFits(await response.text());
    render();
  } catch {
    const warning = node("p", "loading", "The fitted results could not be loaded.");
    const retry = node("button", "", "Try again");
    retry.type = "button";
    retry.onclick = loadFits;
    warning.append(retry);
    plot.replaceChildren(warning);
  }
}
for (let i = 1; i <= 12; i++) byId("replicate").add(new Option("Dataset " + i, String(i)));
form.addEventListener("change", render);
byId("extra").addEventListener("change", render);
for (const mode of ["limits", "error"]) {
  byId("view-" + mode).addEventListener("click", () => {
    display = mode;
    for (const name of ["limits", "error"]) byId("view-" + name).setAttribute("aria-pressed", String(name === mode));
    render();
  });
}
byId("copy").addEventListener("click", async () => {
  try {
    await navigator.clipboard.writeText(byId("install-code").textContent);
    byId("copy-status").textContent = "R code copied.";
    byId("copy").textContent = "Copied";
    window.setTimeout(() => { byId("copy").textContent = "Copy code"; }, 1800);
  } catch {
    const range = document.createRange();
    range.selectNodeContents(byId("install-code"));
    const selection = window.getSelection();
    selection.removeAllRanges();
    selection.addRange(range);
    byId("copy-status").textContent = "Code selected. Use your device’s copy command.";
  }
});
render();
loadFits();
