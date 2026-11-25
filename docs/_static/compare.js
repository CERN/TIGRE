// Landing page: the comparison slider and its tabs.
document.addEventListener("DOMContentLoaded", () => {
  // Move the divider.
  for (const frame of document.querySelectorAll(".tigre-compare-frame")) {
    const range = frame.querySelector(".tigre-compare-range");
    if (!range) continue;
    const update = () => frame.style.setProperty("--pos", `${range.value}%`);
    range.addEventListener("input", update);
    update();
  }
  // Tabs: click, or the arrow, Home and End keys.
  for (const list of document.querySelectorAll(".tigre-compare-tabs")) {
    const tabs = [...list.querySelectorAll('[role="tab"]')];
    const select = (tab, focus) => {
      for (const t of tabs) {
        const on = t === tab;
        t.setAttribute("aria-selected", on);
        t.tabIndex = on ? 0 : -1;
        document.getElementById(t.getAttribute("aria-controls")).hidden = !on;
      }
      if (focus) tab.focus();
    };
    tabs.forEach((tab, i) => {
      tab.addEventListener("click", () => select(tab, false));
      tab.addEventListener("keydown", (e) => {
        const next = { ArrowRight: i + 1, ArrowLeft: i - 1, Home: 0, End: tabs.length - 1 }[e.key];
        if (next === undefined) return;
        e.preventDefault();
        select(tabs[(next + tabs.length) % tabs.length], true);
      });
    });
  }
});
