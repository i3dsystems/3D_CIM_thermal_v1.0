﻿var toc;
(function (toc) {
    var showDetails = false;

    function createOutline(outline) {
        var ol = document.createElement("ol");
        ol.className = "toc";
        for (var i = 0; i < outline.length; i++) {
            ol.appendChild(createSection(outline[i]));
        }
        return ol;
    }
    toc.createOutline = createOutline;

    function createSection(section) {
        var li = document.createElement("li");

        var title = document.createElement("a");
        title.className = "toc_sec_title";
        li.appendChild(title);

        if (section.heading === null) {
            switch (section.associatedNodes[0].nodeName.toLowerCase()) {
                case "blockquote":
                    title.textContent = "Quoted content";
                    break;
                case "body":
                    title.textContent = "Document";
                    break;
                case "details":
                    title.textContent = "Widget";
                    break;
                case "fieldset":
                    title.textContent = "Form controls";
                    break;
                case "figure":
                    title.textContent = "Figure";
                    break;
                case "td":
                    title.textContent = "Data cell";
                    break;
                case "article":
                    title.textContent = "Article";
                    break;
                case "aside":
                    title.textContent = "Aside";
                    break;
                case "nav":
                    title.textContent = "Navigation";
                    break;
                case "section":
                    title.textContent = "Section";
                    break;
            }
            title.className += "toc_no_title";
        } else {
            title.textContent = section.heading.text;
        }

        var node = section.associatedNodes[0];
        if ((node.sectionType !== 1 && node.sectionType !== 2) ||
            (node.nodeName.toLowerCase() === "body")) {
            node = section.heading;
        }
        title.href = "#" + node.id;

        title.addEventListener("click", function (event) {
            event.preventDefault();
            node.scrollIntoView();
        }, false);

        li.appendChild(createOutline(section.childSections));
        return li;
    }

    function Section() {
        this.parentSection = null;
        this.childSections = [];
        this.firstChild = null;
        this.lastChild = null;
        this.appendChild = function (section) {
            section.parentSection = this;
            this.childSections.push(section);
            if (this.firstChild === null) {
                this.firstChild = section;
            }
            this.lastChild = section;
        };

        this.heading = null;

        this.associatedNodes = [];
    }
    toc.Section = Section;

    function HTMLOutline(root) {
        var currentOutlinee = null;

        var currentSection = null;

        var stack = { "lastIndex": -1, "isEmpty": null, "push": null, "pop": null, "top": null };
        stack.isEmpty = function () {
            return stack.lastIndex === -1;
        };
        stack.push = function (e) {
            stack[++stack.lastIndex] = e;
            stack.top = e;
        };
        stack.pop = function () {
            var e = stack.top;
            delete stack[stack.lastIndex--];
            stack.top = stack[stack.lastIndex];
            return e;
        };

        function enter(node) {
            if (isElement(node)) {
                if (!stack.isEmpty() && (isHeadingElement(stack.top) || isHidden(stack.top))) {
                } else if (isHidden(node)) {
                    stack.push(node);
                } else if (isSectioningContentElement(node) || isSectioningRootElement(node)) {
                    if (currentOutlinee !== null) {
                        stack.push(currentOutlinee);
                    }
                    currentOutlinee = node;
                    currentSection = new Section();
                    associateNodeWithSection(currentOutlinee, currentSection);
                    currentOutlinee.appendSection(currentSection);
                } else if (currentOutlinee === null) {
                } else if (isHeadingElement(node)) {
                    if (currentSection.heading === null)
                        currentSection.heading = node;
else if (currentOutlinee.lastSection.heading === null || node.rank >= currentOutlinee.lastSection.heading.rank) {
                        currentSection = new Section();
                        currentSection.heading = node;
                        currentOutlinee.appendSection(currentSection);
                    } else {
                        var candidateSection = currentSection;
                        do {
                            if (node.rank < candidateSection.heading.rank) {
                                currentSection = new Section();
                                currentSection.heading = node;
                                candidateSection.appendChild(currentSection);
                                break;
                            }
                            var newCandidate = candidateSection.parentSection;
                            candidateSection = newCandidate;
                        } while(true);
                    }
                    stack.push(node);
                }
            }
        }

        function exit(node) {
            if (isElement(node)) {
                if (!stack.isEmpty() && node === stack.top) {
                    stack.pop();
                } else if (!stack.isEmpty() && (isHeadingElement(stack.top) || isHidden(stack.top))) {
                } else if (!stack.isEmpty() && isSectioningContentElement(node)) {
                    currentOutlinee = stack.pop();
                    currentSection = currentOutlinee.lastSection;
                    for (var i = 0; i < node.sectionList.length; i++) {
                        currentSection.appendChild(node.sectionList[i]);
                    }
                } else if (!stack.isEmpty() && isSectioningRootElement(node)) {
                    currentOutlinee = stack.pop();
                    currentSection = currentOutlinee.lastSection;
                    while (currentSection.childSections.length > 0) {
                        currentSection = currentSection.lastChild;
                    }
                } else if (isSectioningContentElement(node) || isSectioningRootElement(node)) {
                    currentOutlinee = null;
                    currentSection = null;
                }
            }
            if (node.associatedSection === null && currentSection !== null) {
                associateNodeWithSection(node, currentSection);
            }
        }

        function associateNodeWithSection(node, section) {
            section.associatedNodes.push(node);
            node.associatedSection = section;
        }

        function isElement(node) {
            return node.nodeType === 1;
        }

        function isHidden(node) {
            return node.hidden;
        }

        function isSectioningContentElement(node) {
            return node.sectionType === 1;
        }

        function isSectioningRootElement(node) {
            return node.sectionType === 2;
        }

        function isHeadingElement(node) {
            return node.rank !== undefined;
        }

        function extend(node) {
            if (node.nodeType === 1) {
                switch (node.nodeName.toLowerCase()) {
                    case "blockquote":
                    case "body":
                    case "details":
                    case "dialog":
                    case "fieldset":
                    case "figure":
                    case "td":
                        extendSectioningRootElement(node);
                        break;
                    case "article":
                    case "aside":
                    case "nav":
                    case "section":
                        extendSectioningContentElement(node);
                        break;
                    case "h1":
                    case "h2":
                    case "h3":
                    case "h4":
                    case "h5":
                    case "h6":
                        extendHeadingElement(node);
                        break;
                    case "hgroup":
                        extendHeadingGroupElement(node);
                        break;
                    default:
                        extendNode(node);
                }
            } else
                extendNode(node);
        }

        function extendNode(node) {
            node.associatedSection = null;
        }

        function extendSectioningElement(node) {
            extendNode(node);
            node.sectionList = [];
            node.firstSection = null;
            node.lastSection = null;

            node.appendSection = function (section) {
                this.sectionList.push(section);
                if (this.firstSection === null) {
                    this.firstSection = section;
                }
                this.lastSection = section;
            };
        }

        function extendSectioningContentElement(node) {
            extendSectioningElement(node);
            node.sectionType = 1;
        }

        function extendSectioningRootElement(node) {
            extendSectioningElement(node);
            node.sectionType = 2;
        }

        function extendHeadingContentElement(node) {
            extendNode(node);
            Object.defineProperty(node, "depth", {
                "get": function () {
                    var section = node.associatedSection;
                    var depth = 1;
                    if (section !== null) {
                        while (section = section.parentSection)
                            ++depth;
                    }
                    return depth;
                },
                "configurable": true,
                "enumerable": true
            });
        }

        function extendHeadingElement(node) {
            extendHeadingContentElement(node);
            node.rank = -parseInt(node.nodeName.charAt(1));
            node.text = node.textContent;
        }

        function extendHeadingGroupElement(node) {
            extendHeadingContentElement(node);

            for (var i = 1; i <= 6; i++) {
                var h = node.getElementsByTagName("h" + i);
                if (h.length > 0) {
                    node.rank = -i;
                    node.text = h[0].textContent;
                    break;
                }
            }

            if (node.rank === undefined) {
                node.rank = -1;
                node.text = "";
            }
        }

        var node = root;
        start:
        while (node) {
            extend(node);
            enter(node);
            if (node.firstChild) {
                node = node.firstChild;
                continue start;
            }
            while (node) {
                exit(node);
                if (node === root) {
                    break start;
                }
                if (node.nextSibling) {
                    node = node.nextSibling;
                    continue start;
                }
                node = node.parentNode;
            }
        }
    }
    toc.HTMLOutline = HTMLOutline;
})(toc || (toc = {}));
