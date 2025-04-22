function run() {
    if (!check_path())
	return;

    $('#status').html("<div class=message><div id=msg_logo><img src=/img/gears.gif></div><div id=msg_label>The code is compiling.</div></div>");
    $('#runButton').attr('disabled','disabled');
    $('#cancel').attr('disabled','disabled');
    $('#update').attr('disabled','disabled');

    var url = location.pathname.replace(/_edit\//,"_run/");
    $.post(
	url,
	$("#editedText").serialize(),
	function(data) {
	    $('#messages').html(data);
	    updatePreviewPane();
	    $('#runButton').removeAttr('disabled');
	    $('#cancel').removeAttr('disabled');
	    $('#update').removeAttr('disabled');
	},
	"html");
};

$(document).ready(function(){
    $("#runButton").show();
});
