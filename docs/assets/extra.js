/* Keep the search results in the language of the page being read.
 *
 * mkdocs-static-i18n builds one search index for the whole site: its
 * `extend_search_entries` stacks every locale's entries together, with no
 * option to separate them, and Material fetches that one file relative to the
 * site root from both locales. So an English reader searching `gff` got every
 * matching page twice, once in English and once in Spanish, and a Spanish
 * reader got the mirror image. On a reference site most queries are made of
 * language-neutral words (gff, vcf, bam, indel, a column name, a flag), so this
 * was most queries.
 *
 * There is no server-side fix inside the plugin, so the results that do not
 * belong to this page's language are dropped as they are rendered. The filter
 * reads the language from the document and the site root from Material's own
 * config, so a third locale needs nothing here beyond its prefix.
 *
 * It is keyed to Material's result classes. If a future release renames them
 * the observer simply never fires and the site falls back to showing both
 * languages, which is today's behaviour: it fails open, not broken.
 */
(function () {
  var config = document.getElementById("__config");
  if (!config) return;

  var root;
  try {
    root = new URL(JSON.parse(config.textContent).base, location.href).pathname;
  } catch (error) {
    return;
  }
  if (!root.endsWith("/")) root += "/";

  // Every locale that is built as a subdirectory. The default locale has no
  // prefix, so anything that is not one of these belongs to it.
  var PREFIXES = ["es"];
  var language = document.documentElement.lang || "en";

  var list = document.querySelector(".md-search-result__list");
  var meta = document.querySelector(".md-search-result__meta");
  if (!list) return;

  var localeOf = function (href) {
    var segment = new URL(href, location.href).pathname.slice(root.length).split("/")[0];
    return PREFIXES.indexOf(segment) >= 0 ? segment : "en";
  };

  var pruning = false;
  new MutationObserver(function () {
    if (pruning) return;
    pruning = true;
    var kept = 0;
    list.querySelectorAll("li.md-search-result__item").forEach(function (item) {
      var link = item.querySelector("a.md-search-result__link");
      if (!link) return;
      if (localeOf(link.href) === language) kept++;
      else item.remove();
    });
    // The counter now says how many results are shown rather than how many
    // index entries matched. Both locales' strings lead with the number.
    if (meta && /\d/.test(meta.textContent)) {
      meta.textContent = meta.textContent.replace(/\d+/, String(kept));
    }
    pruning = false;
  }).observe(list, { childList: true });
})();
