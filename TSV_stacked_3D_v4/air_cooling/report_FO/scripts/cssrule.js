var cssrule;
(function (cssrule) {

    function createStyleSheet() {
        var style = document.createElement("style");

        // Add a media (and/or media query) here if you'd like!
        // style.setAttribute("media", "screen")
        // style.setAttribute("media", "@media only screen and (max-width : 1024px)")

        // WebKit hack :(
        style.appendChild(document.createTextNode(""));

        // Add the <style> element to the page
        document.head.appendChild(style);

        return style.sheet;
    }

    function addCSSRule(selector, rules, index) {
        if (document.styleSheets) {
            var sheet = createStyleSheet();
            if (sheet) {
                if (sheet.insertRule) {
                    return sheet.insertRule(selector + " {" + rules + "}", index);
                }
                else {
                    return sheet.addRule(selector, rules, index);
                }
            }
        }
        return null;
    }
    cssrule.addCSSRule = addCSSRule;

})(cssrule || (cssrule = {}));