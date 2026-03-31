(function(window, $) {
    function shouldEnableGeneAutocomplete(term) {
        const value = String(term || "").trim();

        if (value.length < 3) {
            return false;
        }

        return !/[\s:;]/.test(value);
    }

    function initGeneAutocomplete(inputSelector, availableTags, options) {
        const settings = options || {};
        const $input = $(inputSelector);

        if ($input.length === 0 || typeof $input.autocomplete !== "function") {
            return;
        }

        const tags = Array.isArray(availableTags) ? availableTags : [];

        $input.autocomplete({
            minLength: 3,
            delay: 500,
            source: function(request, response) {
                const term = request.term || "";

                if (!shouldEnableGeneAutocomplete(term)) {
                    response([]);
                    return;
                }

                response($.ui.autocomplete.filter(tags, term));
            }
        });

        $input.on("input", function() {
            const currentValue = $(this).val() || "";
            if (!shouldEnableGeneAutocomplete(currentValue)) {
                $input.autocomplete("close");
            }
        });

        if (settings.closeOnBlur !== false) {
            $input.on("blur", function() {
                $input.autocomplete("close");
            });
        }
    }

    window.METADOME_UI = window.METADOME_UI || {};
    window.METADOME_UI.shouldEnableGeneAutocomplete = shouldEnableGeneAutocomplete;
    window.METADOME_UI.initGeneAutocomplete = initGeneAutocomplete;
})(window, jQuery);