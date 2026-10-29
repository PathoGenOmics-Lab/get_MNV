/* The home page's two moving parts: Figure 1, which steps through one real
 * codon of the bundled sample, and the codon playground. Plain script, no
 * build step, works from file://. Without it the page still reads: the figure
 * shows its last step and the playground its first example, worked out. */
(function () {
  var es = document.documentElement.lang === "es";

  // ── Figure 1 ──

  var fig = document.querySelector("[data-fig]");
  if (fig) {
    var dots = Array.prototype.slice.call(fig.querySelectorAll("[data-go]"));
    var count = fig.querySelector("[data-fig-count]");
    var title = fig.querySelector("[data-fig-title]");
    var play = fig.querySelector("[data-fig-play]");
    var steps = dots.length;
    var timer = null;
    var still = window.matchMedia && matchMedia("(prefers-reduced-motion: reduce)").matches;

    var titleOf = function (n) {
      var label = dots[n - 1].getAttribute("aria-label");
      var text = label.slice(label.indexOf(":") + 1).trim();
      return text.charAt(0).toUpperCase() + text.slice(1);
    };
    var show = function (n) {
      fig.dataset.step = String(n);
      dots.forEach(function (d, i) { d.setAttribute("aria-pressed", i + 1 === n ? "true" : "false"); });
      count.textContent = n + (es ? " de " : " of ") + steps;
      title.textContent = titleOf(n);
    };
    var stop = function () {
      clearTimeout(timer);
      timer = null;
      play.textContent = es ? "Reproducir" : "Play";
    };
    var run = function () {
      play.textContent = es ? "Pausa" : "Pause";
      var at = +fig.dataset.step;
      timer = setTimeout(function () {
        show(at >= steps ? 1 : at + 1);
        run();
      }, at >= steps ? 7000 : 4200);
    };

    dots.forEach(function (d, i) {
      d.addEventListener("click", function () { stop(); show(i + 1); });
      d.addEventListener("keydown", function (event) {
        if (event.key !== "ArrowRight" && event.key !== "ArrowLeft") return;
        event.preventDefault();
        var next = (i + (event.key === "ArrowRight" ? 1 : steps - 1)) % steps;
        stop();
        show(next + 1);
        dots[next].focus();
      });
    });
    play.hidden = false;
    play.addEventListener("click", function () { if (timer) stop(); else run(); });

    // A reader who asked for less motion gets the finished figure and the
    // steps to walk through by hand; everyone else gets it told from step 1.
    if (still) {
      show(steps);
      stop();
    } else {
      show(1);
      run();
    }

    // The same position in the VCF and in the gene, lit together.
    fig.querySelectorAll("[data-pos]").forEach(function (el) {
      var pos = el.getAttribute("data-pos");
      var twins = fig.querySelectorAll('[data-pos="' + pos + '"]');
      el.addEventListener("mouseenter", function () { twins.forEach(function (t) { t.classList.add("lit"); }); });
      el.addEventListener("mouseleave", function () { twins.forEach(function (t) { t.classList.remove("lit"); }); });
    });
  }

  // ── The codon playground ──

  var playground = document.querySelector("[data-play]");
  if (!playground) return;

  var BASES = "TCAG";
  // The standard code, which table 11 shares for every codon inside a gene.
  var CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  var NAME = {
    A: "Ala", R: "Arg", N: "Asn", D: "Asp", C: "Cys", Q: "Gln", E: "Glu", G: "Gly", H: "His", I: "Ile",
    L: "Leu", K: "Lys", M: "Met", F: "Phe", P: "Pro", S: "Ser", T: "Thr", W: "Trp", Y: "Tyr", V: "Val",
    "*": "Stop",
  };
  var T = es ? {
    head: ["Leído como", "Codón", "Aminoácido"],
    reference: "Referencia",
    alone: "Posición {i} sola, {r}>{a}",
    single: "Posición {i}, {r}>{a}",
    whole: ["", "", "Las dos, como lo lee get_MNV", "Las tres, como lo lee get_MNV"],
    kinds: { syn: "sinónima", mis: "missense", gain: "stop ganado", lost: "stop perdido" },
    position: "Posición {i}",
    group: "Base en la posición {i}",
    refTag: "ref.",
    refInput: "Codón de referencia",
    and: " y ",
    none: "Elige una base para cambiarla.",
    one: "Con un solo cambio en el codón, leer la base y leer el codón coinciden.",
    differs: ["", "", "Base a base, este codón da {list}. Leído entero da {aa}, que no nombra ninguna de las dos.",
      "Base a base, este codón da {list}. Leído entero da {aa}, que no nombra ninguna de ellas."],
    masked: "Por sí sola, la posición {i} da un stop. Junto {others}, el codón da {aa} y la proteína sigue.",
    others: ["", "", "al otro cambio", "a los otros cambios"],
    gained: "Cada cambio es silencioso por separado. Juntos dan {aa}.",
    stopTogether: ["", "", "Ningún cambio da un stop por sí solo; juntos, sí.", "Ningún cambio da un stop por sí solo; juntos, sí."],
    same: "Leído entero da {aa}, como la posición {i} por sí sola.",
  } : {
    head: ["Read as", "Codon", "Amino acid"],
    reference: "Reference",
    alone: "Position {i} alone, {r}>{a}",
    single: "Position {i}, {r}>{a}",
    whole: ["", "", "Both, as get_MNV reads it", "All three, as get_MNV reads it"],
    kinds: { syn: "synonymous", mis: "missense", gain: "stop gained", lost: "stop lost" },
    position: "Position {i}",
    group: "Base at position {i}",
    refTag: "ref",
    refInput: "Reference codon",
    and: " and ",
    none: "Pick a base to change it.",
    one: "With one change in the codon, reading the base and reading the codon agree.",
    differs: ["", "", "One base at a time this codon gives {list}. Read whole, it gives {aa}, which neither of them names.",
      "One base at a time this codon gives {list}. Read whole, it gives {aa}, which none of them names."],
    masked: "On its own, position {i} makes a stop. Together with {others} the codon makes {aa}, and the protein goes on.",
    others: ["", "", "the other change", "the other changes"],
    gained: "Each change is silent on its own. Together they make {aa}.",
    stopTogether: ["", "", "Neither change makes a stop on its own; together they do.", "No single change makes a stop; together they do."],
    same: "Read whole, it gives {aa}, as position {i} does on its own.",
  };
  var fill = function (text, values) {
    return text.replace(/\{(\w+)\}/g, function (_, key) { return values[key]; });
  };
  var aaOf = function (codon) {
    var i = BASES.indexOf(codon[0]) * 16 + BASES.indexOf(codon[1]) * 4 + BASES.indexOf(codon[2]);
    return CODE[i];
  };
  var kindOf = function (ref, alt) {
    if (ref === alt) return "syn";
    if (alt === "*") return "gain";
    if (ref === "*") return "lost";
    return "mis";
  };
  var listOf = function (names) {
    return names.length < 2 ? names.join("") : names.slice(0, -1).join(", ") + T.and + names[names.length - 1];
  };

  var ref = ["C", "G", "T"];
  var cur = ["T", "C", "T"];
  var editor = playground.querySelector("[data-editor]");
  var out = playground.querySelector("[data-out]");
  var presets = Array.prototype.slice.call(playground.querySelectorAll("[data-preset]"));

  // Everything below is built as nodes with text in them, never as markup: the
  // bases come from an input box and from attributes, and none of it is ever
  // read as HTML, whatever it holds.
  var el = function (tag, attrs, children) {
    var node = document.createElement(tag);
    Object.keys(attrs || {}).forEach(function (key) { node.setAttribute(key, attrs[key]); });
    (children || []).forEach(function (child) {
      node.appendChild(typeof child === "string" ? document.createTextNode(child) : child);
    });
    return node;
  };

  // The editor: three positions, each with the four bases to choose from.
  editor.textContent = "";
  var picks = [];
  var tiles = [];
  [0, 1, 2].forEach(function (p) {
    var col = document.createElement("div");
    col.className = "pos";
    var label = document.createElement("p");
    label.className = "pos-label";
    label.textContent = fill(T.position, { i: p + 1 });
    var tile = document.createElement("div");
    tile.className = "nt big";
    tile.setAttribute("aria-hidden", "true");
    var group = document.createElement("div");
    group.className = "pick";
    group.setAttribute("role", "radiogroup");
    group.setAttribute("aria-label", fill(T.group, { i: p + 1 }));
    var buttons = "ACGT".split("").map(function (b) {
      var button = document.createElement("button");
      button.type = "button";
      button.className = "nt";
      button.dataset.b = b;
      button.setAttribute("role", "radio");
      button.textContent = b;
      button.addEventListener("click", function () { cur[p] = b; choose(null); });
      button.addEventListener("keydown", function (event) {
        var at = "ACGT".indexOf(b);
        var step = { ArrowRight: 1, ArrowDown: 1, ArrowLeft: 3, ArrowUp: 3 }[event.key];
        if (!step) return;
        event.preventDefault();
        cur[p] = "ACGT"[(at + step) % 4];
        choose(null);
        picks[p][(at + step) % 4].focus();
      });
      group.appendChild(button);
      return button;
    });
    picks.push(buttons);
    tiles.push(tile);
    col.appendChild(label);
    col.appendChild(tile);
    col.appendChild(group);
    editor.appendChild(col);
  });
  var refInput = el("input", {
    type: "text", inputmode: "text", maxlength: "3", spellcheck: "false",
    autocomplete: "off", autocapitalize: "characters",
  });
  var refLabel = el("label", { class: "play-ref" }, [el("span", null, [T.refInput]), refInput]);
  refInput.addEventListener("input", function () {
    var value = refInput.value.toUpperCase().replace(/U/g, "T");
    refInput.setAttribute("aria-invalid", /^[ACGT]{3}$/.test(value) ? "false" : "true");
    if (!/^[ACGT]{3}$/.test(value)) return;
    ref = value.split("");
    cur = ref.slice();
    choose(null);
  });
  editor.appendChild(refLabel);

  function render() {
    var refAa = aaOf(ref.join(""));
    var changed = [0, 1, 2].filter(function (p) { return cur[p] !== ref[p]; });
    var whole = cur.join("");
    var wholeAa = aaOf(whole);

    picks.forEach(function (buttons, p) {
      buttons.forEach(function (button) {
        var on = button.dataset.b === cur[p];
        button.setAttribute("aria-checked", on ? "true" : "false");
        button.tabIndex = on ? 0 : -1;
        button.classList.toggle("is-ref", button.dataset.b === ref[p]);
        button.title = button.dataset.b === ref[p] ? T.refTag : "";
      });
      tiles[p].dataset.b = cur[p];
      tiles[p].classList.toggle("changed", cur[p] !== ref[p]);
      tiles[p].textContent = cur[p];
      tiles[p].appendChild(el("small", null, [cur[p] !== ref[p] ? ref[p] + ">" + cur[p] : T.refTag]));
    });
    if (document.activeElement !== refInput) refInput.value = ref.join("");

    var codonNode = function (codon) {
      return el("code", null, codon.split("").map(function (b, p) {
        return b !== ref[p] ? el("i", null, [b]) : b;
      }));
    };
    var aaNodes = function (aa) {
      var kind = kindOf(refAa, aa);
      return [NAME[aa] + " ", el("span", { class: "kind kind-" + kind }, [T.kinds[kind]])];
    };
    var row = function (label, codon, aa, whole) {
      return el("tr", whole ? { class: "whole" } : null, [
        el("th", { scope: "row" }, [label]),
        el("td", null, [codonNode(codon)]),
        el("td", null, aa),
      ]);
    };
    var rows = [row(T.reference, ref.join(""), [NAME[refAa]], false)];
    var singles = changed.map(function (p) {
      var codon = ref.slice();
      codon[p] = cur[p];
      codon = codon.join("");
      var aa = aaOf(codon);
      var label = fill(changed.length > 1 ? T.alone : T.single, { i: p + 1, r: ref[p], a: cur[p] });
      rows.push(row(label, codon, aaNodes(aa), changed.length === 1));
      return { p: p, aa: aa };
    });
    if (changed.length > 1) {
      rows.push(row(T.whole[changed.length], whole, aaNodes(wholeAa), true));
    }

    var verdict;
    if (!changed.length) verdict = T.none;
    else if (changed.length === 1) verdict = T.one;
    else {
      var stopAlone = singles.filter(function (s) { return s.aa === "*" && refAa !== "*"; })[0];
      var silent = singles.every(function (s) { return s.aa === refAa; });
      var match = singles.filter(function (s) { return s.aa === wholeAa; })[0];
      if (stopAlone && wholeAa !== "*") verdict = fill(T.masked, { i: stopAlone.p + 1, aa: NAME[wholeAa], others: T.others[changed.length] });
      else if (silent && wholeAa !== refAa) verdict = fill(T.gained, { aa: NAME[wholeAa] });
      else if (wholeAa === "*" && refAa !== "*" && !stopAlone) verdict = T.stopTogether[changed.length];
      else if (match) verdict = fill(T.same, { aa: NAME[wholeAa], i: match.p + 1 });
      else {
        var names = singles.map(function (s) { return NAME[s.aa]; });
        verdict = fill(T.differs[changed.length], { list: listOf(names), aa: NAME[wholeAa] });
      }
    }

    var head = el("thead", null, [el("tr", null, T.head.map(function (h) { return el("th", { scope: "col" }, [h]); }))]);
    out.textContent = "";
    out.appendChild(el("table", { class: "readout" }, [head, el("tbody", null, rows)]));
    out.appendChild(el("p", { class: "verdict" }, [verdict]));
  }

  function choose(preset) {
    presets.forEach(function (b) { b.setAttribute("aria-pressed", b === preset ? "true" : "false"); });
    render();
  }

  presets.forEach(function (button) {
    button.addEventListener("click", function () {
      var parts = button.getAttribute("data-preset").split(" ");
      ref = parts[0].split("");
      cur = parts[1].split("");
      choose(button);
    });
  });

  playground.classList.add("live");
  choose(presets[0]);
})();
