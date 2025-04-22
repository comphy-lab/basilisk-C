function toline(line)
{
    var t = editor.charCoords ({line: line - 1, ch: 0}, "local").top; 
    editor.scrollTo (null, t);
}

function check_path() {
    checked = true;
    $('input[type=checkbox]').each(function () {
	if (!this.checked)
	    checked = false;
    });
    if (!checked)
	$('#messages').html("<li>Please confirm the creation of new directories.</li>");
    return checked;
}

function save() {
    if (!check_path())
	return;
    var url = location.pathname.replace(/_edit\//,"_save/");
    $.post(url, $("#logMsg, #editedText").serialize(),
	   function(data) {		   
	       console.log(data);
	       if (data != "")
		   $('#messages').html(data);
	       else
		   window.location.href = location.href.replace(/\/_edit\//,"/");
	   }
	  );
}

function checkRunning(url, sec) {
    if ($('#status').html().search("running") != -1) {
	setTimeout (function() {
	    $('#status').load(url, function(data) {
		remaining = data.match(/([0-9]+):([0-9]+)/);
		if (remaining)
		    sec = Math.max (100*(parseFloat(remaining[1])*60. +
					 parseFloat(remaining[2])), 2000);
		else
		    sec = 2*sec;
		checkRunning (url, Math.min(sec, 16000));
	    });
	}, sec);
    }
    else if (sec > 1000)
	updatePreviewPane(1);
}

function updatePreviewPane(nostatus) {
    $('#previewButton').attr('disabled','disabled');
    var url = location.pathname.replace(/_edit\//,"_preview/");
    $.post(
	url,
	{"raw" : $("#editedText").attr("value")},
	function(data) {
            $('#previewpane').html(data);
	    hide_plots();
	    render_math();
	    if (!nostatus && location.pathname.match ("[.][cm]$")) {
		var url = location.pathname.concat("?status");
		$('#status').load(url, function(data) { checkRunning (url, 1000); });
	    }
	    $('#previewButton').removeAttr('disabled');
	},
	"html");
};

$(document).ready(function() {
    $("#previewButton").show();
    var url = location.pathname.replace(/_edit\//,"") + "?raw";
    $.post (url, function(data) {
	if (data.match("^<h1>404 Not Found</h1>"))
	    $("#editedText").attr("value", "");
	else
	    $("#editedText").attr("value", data);
	setup_codemirror();
	$("#editedText").focus();
	updatePreviewPane();
    });
});

$(window).unload(function(){
    var request = new XMLHttpRequest();
    request.open ("POST", location.pathname, false);
    request.setRequestHeader ("content-type", "application/x-www-form-urlencoded");
    request.send ("cancel=Discard");
});
