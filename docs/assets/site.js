/* The documentation site's behaviour: the theme switch, the search box in the
 * header, copy buttons on code, the menu on small screens, and the table of
 * contents following the reader. Plain script, no build step, works from
 * file://. The home page's figure and playground are in home.js. */
(function () {
  var root = document.documentElement;

  // Theme: the reader's choice is kept; until they make one, dark.
  var toggle = document.querySelector("[data-theme-toggle]");
  if (toggle) {
    toggle.addEventListener("click", function () {
      var next = root.dataset.theme === "dark" ? "light" : "dark";
      root.dataset.theme = next;
      try { localStorage.setItem("get_mnv.theme", next); } catch (e) {}
    });
  }

  // Copy buttons on code blocks.
  var copyLabel = root.lang === "es" ? "Copiar" : "Copy";
  var copiedLabel = root.lang === "es" ? "Copiado" : "Copied";
  document.querySelectorAll(".code").forEach(function (block) {
    var button = document.createElement("button");
    button.type = "button";
    button.className = "copy";
    button.textContent = copyLabel;
    button.addEventListener("click", function () {
      var text = block.querySelector("pre").innerText;
      var done = function () {
        button.textContent = copiedLabel;
        setTimeout(function () { button.textContent = copyLabel; }, 1500);
      };
      if (navigator.clipboard) navigator.clipboard.writeText(text).then(done, function () {});
    });
    block.appendChild(button);
  });

  // The page menu starts folded on a phone, where it would push the page down.
  var menu = document.querySelector(".menu");
  if (menu && window.matchMedia("(max-width: 820px)").matches) menu.removeAttribute("open");

  // Search over the index built with the pages, in the language of this page.
  var input = document.querySelector("[data-search]");
  var box = document.querySelector("[data-results]");
  var index = (window.GET_MNV_DOCS_INDEX || []).filter(function (e) { return e.lang === root.lang; });
  var base = root.dataset.root || "";
  var esc = function (s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c];
    });
  };
  var picked = -1;
  function search(query) {
    var words = query.toLowerCase().split(/\s+/).filter(Boolean);
    if (!words.length) return [];
    var hits = [];
    index.forEach(function (e) {
      var head = (e.page + " " + e.section).toLowerCase();
      var body = e.text.toLowerCase();
      var score = 0;
      for (var i = 0; i < words.length; i++) {
        var inHead = head.indexOf(words[i]) >= 0, inBody = body.indexOf(words[i]) >= 0;
        if (!inHead && !inBody) return;
        score += inHead ? 3 : 1;
      }
      hits.push({ e: e, score: score });
    });
    hits.sort(function (a, b) { return b.score - a.score; });
    return hits.slice(0, 12);
  }
  function snippet(text, word) {
    var at = text.toLowerCase().indexOf(word);
    if (at < 0) return text.slice(0, 120);
    var start = Math.max(0, at - 50);
    return (start ? "…" : "") + text.slice(start, start + 140) + "…";
  }
  function show() {
    var query = input.value.trim();
    picked = -1;
    if (!query) { box.hidden = true; box.innerHTML = ""; return; }
    var hits = search(query);
    var first = query.toLowerCase().split(/\s+/)[0];
    box.innerHTML = hits.length
      ? hits.map(function (h) {
          var title = h.e.section ? h.e.page + " · " + h.e.section : h.e.page;
          return '<a href="' + esc(base + h.e.url) + '"><b>' + esc(title) + "</b><span>" + esc(snippet(h.e.text, first)) + "</span></a>";
        }).join("")
      : '<p class="none">' + (root.lang === "es" ? "Sin resultados" : "No results") + "</p>";
    box.hidden = false;
  }
  if (input && box) {
    input.addEventListener("input", show);
    input.addEventListener("keydown", function (event) {
      var links = box.querySelectorAll("a");
      if (event.key === "ArrowDown" || event.key === "ArrowUp") {
        if (!links.length) return;
        event.preventDefault();
        // From nothing picked, down goes to the first result and up to the last.
        if (picked < 0) picked = event.key === "ArrowDown" ? 0 : links.length - 1;
        else picked = (picked + (event.key === "ArrowDown" ? 1 : -1) + links.length) % links.length;
        links.forEach(function (a, i) { a.classList.toggle("on", i === picked); });
        links[picked].scrollIntoView({ block: "nearest" });
      } else if (event.key === "Enter" && links.length) {
        location.href = links[Math.max(0, picked)].href;
      } else if (event.key === "Escape") {
        input.value = ""; show();
        if (field) field.classList.remove("open");
        input.blur();
      }
    });
    document.addEventListener("click", function (event) {
      if (!box.contains(event.target) && event.target !== input) box.hidden = true;
      if (field && !field.contains(event.target) && !input.value) field.classList.remove("open");
    });
    // On a phone the box opens from its button, over the page; anywhere, "/"
    // puts the cursor in it, unless the reader is already typing somewhere.
    var field = document.querySelector(".search");
    var opener = document.querySelector("[data-search-open]");
    if (opener && field) {
      opener.addEventListener("click", function () { field.classList.add("open"); input.focus(); });
    }
    document.addEventListener("keydown", function (event) {
      if (event.key !== "/" || event.metaKey || event.ctrlKey || event.altKey) return;
      var target = event.target;
      if (target && (/^(INPUT|TEXTAREA|SELECT)$/.test(target.tagName) || target.isContentEditable)) return;
      event.preventDefault();
      if (field) field.classList.add("open");
      input.focus();
    });
  }

  // The table of contents marks the section being read.
  var tocLinks = Array.prototype.slice.call(document.querySelectorAll(".toc a"));
  if (tocLinks.length && "IntersectionObserver" in window) {
    var byId = {};
    tocLinks.forEach(function (a) { byId[a.getAttribute("href").slice(1)] = a; });
    var observer = new IntersectionObserver(function (seen) {
      seen.forEach(function (entry) {
        if (!entry.isIntersecting) return;
        tocLinks.forEach(function (a) { a.classList.remove("on"); });
        var a = byId[entry.target.id];
        if (a) a.classList.add("on");
      });
    }, { rootMargin: "-60px 0px -70% 0px" });
    Object.keys(byId).forEach(function (id) {
      var heading = document.getElementById(id);
      if (heading) observer.observe(heading);
    });
  }
})();
