// Landing page: the comparison slider and its tabs.
document.addEventListener("DOMContentLoaded", () => {
  // Move the divider.
  for (const frame of document.querySelectorAll(".tigre-compare-frame")) {
    const range = frame.querySelector(".tigre-compare-range");
    if (!range) continue;
    const left = frame.querySelector(".tigre-compare-left");
    const right = frame.querySelector(".tigre-compare-right");
    const update = () => {
      frame.style.setProperty("--pos", `${range.value}%`);
      // The zoom box is anchored at the right edge: cut it where the divider is.
      frame.style.setProperty("--cut", `${(frame.clientWidth * (100 - range.value)) / 100}px`);
      // Hide a name when its side of the divider is too narrow for it.
      const x = (frame.clientWidth * range.value) / 100;
      if (left) left.classList.toggle("is-hidden", x < left.offsetWidth + 16);
      if (right) right.classList.toggle("is-hidden", frame.clientWidth - x < right.offsetWidth + 16);
    };
    range.addEventListener("input", update);
    window.addEventListener("resize", update);
    // A magnified box grows while the divider is held or dragged, and shrinks
    // back a moment after it is let go.
    const zoom = frame.querySelector(".tigre-compare-zoom");
    if (zoom) {
      let timer;
      const grow = () => { clearTimeout(timer); frame.classList.add("is-zoomed"); };
      const shrink = () => { clearTimeout(timer); timer = setTimeout(() => frame.classList.remove("is-zoomed"), 1500); };
      range.addEventListener("pointerdown", grow);
      range.addEventListener("input", grow);
      range.addEventListener("pointerup", shrink);
      range.addEventListener("pointercancel", shrink);
      range.addEventListener("change", shrink);
      range.addEventListener("blur", shrink);
    }
    frame.tigreUpdate = update;
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
        const panel = document.getElementById(t.getAttribute("aria-controls"));
        panel.hidden = !on;
        // A hidden panel has no size, so measure its names once it is shown.
        if (on) panel.querySelectorAll(".tigre-compare-frame").forEach((f) => f.tigreUpdate?.());
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
