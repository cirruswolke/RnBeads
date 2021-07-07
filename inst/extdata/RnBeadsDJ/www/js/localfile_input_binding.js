(function() {
/**
 * Shiny Registration
 * thanks to https://github.com/wleepang/shiny-directory-input
 */

var localFileInputBinding = new Shiny.InputBinding();
$.extend(localFileInputBinding, {
  find: function(scope) {
    return( $(scope).find(".localfile-input") );
  },
  initialize: function(el) {
    // called when document is ready using initial values defined in ui.R
    // documented in input_binding.js but not in docs (articles)
  },
  getId: function(el) {
    return($(el).attr('id'));
  },
  getValue: function(el) {
    return($(el).data('val') || 0);
  },
  setValue: function(el, value) {
    $(el).data('val', value);
  },
  receiveMessage: function(el, data) {
    // This is used for receiving messages that tell the input object to do
    // things, such as setting values (including min, max, and others).
    // documented in input_binding.js but not in docs (articles)
    var $widget = $(el).parentsUntil('.localfile-input-container').parent();
    var $path = $widget.find('input.localfile-input-chosen-file');

    console.log('message received: ' + data.chosen_file);

    if (data.chosen_file) {
      $path.val(data.chosen_file);
      $path.trigger('change');
    }
  },
  subscribe: function(el, callback) {
    $(el).on("click.localFileInputBinding", function(e) {
      var $el = $(this);
      var val = $el.data('val') || 0;
      $el.data('val', val + 1);

      console.log('in subscribe: click');
      callback();
    });
  },
  unsubscribe: function(el) {
    $(el).off(".localFileInputBinding");
  }
});

Shiny.inputBindings
  .register(localFileInputBinding, "oddhypothesis.localFileInputBinding");


})();
