// "use strict";

// Constants
var MAX_UPLOAD_FILE_SIZE = 100*1024*1024; // 100 MB
var UPLOAD_URL = "";
var NEXT_URL   = "/results/";

// File to handle when the submit button is finally clicked.
var dropped_file  = [];

$(document).ready(function() {

  // Set up the drag/drop zone.
  setup_dropzone();

  // Set up the handler for the file input box.
  // Propagate changes to the dropzone text fields
  var $dropzone = $('#drag_n_drop');
  var $input = $('#structure');
  $input.on('change', function() {
    if ($(this).val()) {
      $dropzone.find('.file_upload_text').css('display', 'none');
      $dropzone.find('.drag_n_drop_text').css('display', 'none');
      $dropzone.find('.uploaded_file_name').css('display', 'inline');
      $dropzone.find('.uploaded_file_name').text($(this).val().split('\\').pop());

      // dropped_file = this.files;
      file_list = this.files;
      for (var i = 0, ie = file_list.length; i < ie; i++) {
          dropped_file.push(file_list[i]);
      };

    } else {
      // Nothing selected
      $dropzone.find('.file_upload_text').css('display', 'inline');
      $dropzone.find('.drag_n_drop_text').css('display', 'inline');
      $dropzone.find('.uploaded_file_name').css('display', 'none');
      $dropzone.find('.uploaded_file_name').text('');

      dropped_file = null;
    };
  });

  // Handle the submit button.
  $("#label_submit").on("click", function(e) {
      // If the user has JS disabled, none of this code is running but the
      // file upload input box should still work. In this case they'll
      // just POST to the upload endpoint directly. However, with JS we'll do
      // the POST using ajax and then redirect them ourself when done.
      e.preventDefault();
      do_upload();
  });

  function do_upload() {
      var $progressBar = $("#progress_bar");

      // Gray out the form.
      $("#drag_n_drop :input").attr("disabled", "disabled");

      // Initialize the progress bar.
      $progressBar.css({"width": "0%"});

      // Collect the form data.
      var fd = collectFormData();

      // Attach the file.
      for (var i = 0, ie = dropped_file.length; i < ie; i++) {
          // Collect the other form data.
          fd.append("structure", dropped_file[i]);
      }

      // Inform the back-end that we're doing this over ajax.
      fd.append("__ajax", "true");

      var xhr = $.ajax({
          xhr: function() {
              var xhrobj = $.ajaxSettings.xhr();
              if (dropped_file.length != 0 && xhrobj.upload) {
                  xhrobj.upload.addEventListener("progress", function(event) {
                      var percent = 0;
                      var position = event.loaded || event.position;
                      var total    = event.total;
                      if (event.lengthComputable) {
                          percent = Math.ceil(position / total * 100);
                      };

                      // Set the progress bar.
                      $progressBar.css({"width": percent + "%"});
                  }, false)
              };
              return xhrobj;
          },

          url: UPLOAD_URL,
          method: "POST",
          contentType: false,
          processData: false,
          cache: false,
          data: fd,
          success: function(data) {
              data = JSON.parse(data);

              // How'd it go?
              if (data.status === "ERROR") {
                  // Uh-oh.
                  // console.log(data.msg);
                  $("#drag_n_drop :input").removeAttr("disabled");
                  window.location.reload()
                  return;

              } else {
                  var $dropzone = $('#drag_n_drop');
                  $progressBar.css({"width": "100%"});
                  $dropzone.find('.uploaded_file_name').text('Structure uploaded!');

                  // Ok! Get the Job Name.
                  // Redirect
                  var job_name = data.msg;
                  window.location = NEXT_URL + job_name;
              }
          },
      });
  };

  function collectFormData() {
      // Go through all the form fields and collect their names/values.
      var fd = new FormData();

      $("#drag_n_drop :input").each(function() {
          var $this = $(this);
          var name  = $this.attr("name");
          var type  = $this.attr("type") || "";
          var value = $this.val();

          // No name = no care.
          if (name === undefined) {
              return;
          };

          // Skip the file upload box for now.
          if (type === "file") {
              return;
          };

          // Checkboxes? Only add their value if they're checked.
          if (type === "checkbox" || type === "radio") {
              if (!$this.is(":checked")) {
                  return;
              }
          };

          fd.append(name, value);
      });

      return fd;
  };

  function setup_dropzone() {

    // Partially adapted from
    // github.com/kirsle/flask-multi-upload
    // css-tricks.com/drag-and-drop-file-uploading/
    var $dropzone = $('#drag_n_drop');

    // Show the drag-n-drop field
    $dropzone.on('drag dragstart dragend dragover dragenter dragleave drop', function(e) {
      e.preventDefault();
      e.stopPropagation();
    })

    .on('dragover dragenter', function() {
      $dropzone.removeClass('no_dragover');
      $dropzone.addClass('is_dragover');
    })
    .on('dragleave dragend drop', function() {
      $dropzone.removeClass('is_dragover');
      $dropzone.addClass('no_dragover');
    })

    .on('drop', function(e) {
      // Multiple files can be uploaded, we only want the last.
      file_list = e.originalEvent.dataTransfer.files;
      for (var i = 0, ie = file_list.length; i < ie; i++) {
          dropped_file.push(file_list[i]);
      };

      // Change text
      $dropzone.find('.file_upload_text').css('display', 'none');
      $dropzone.find('.drag_n_drop_text').css('display', 'none');
      $dropzone.find('.uploaded_file_name').css('display', 'inline');
      $dropzone.find('.uploaded_file_name').text(file_list[0].name);
    });

    // If the files are dropped outside of the drop zone, the browser will
    // redirect to show the files in the window. To avoid that we can prevent
    // the 'drop' event on the document.
    function stopDefault(e) {
        e.stopPropagation();
        e.preventDefault();
    };
    $(document).on("dragenter", stopDefault);
    $(document).on("dragover", stopDefault);
    $(document).on("drop", stopDefault);
  };

});
