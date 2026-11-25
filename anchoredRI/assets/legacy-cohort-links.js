// Keep links to sections and controls from the former home page working.
const fragment = window.location.hash.slice(1);
if (fragment && document.getElementById("former-comparison")?.querySelector(
  `[id="${CSS.escape(fragment)}"]`
)) {
  window.location.replace("prepared-cohort.html" + window.location.search + window.location.hash);
}
