// Landing page: move the divider of the FDK/CGLS comparison.
document.addEventListener("DOMContentLoaded", () => {
  for (const frame of document.querySelectorAll(".tigre-compare-frame")) {
    const range = frame.querySelector(".tigre-compare-range");
    if (!range) continue;
    const update = () => frame.style.setProperty("--pos", `${range.value}%`);
    range.addEventListener("input", update);
    update();
  }
});
