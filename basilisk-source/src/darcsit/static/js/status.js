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
	location.reload (true);
}

$(document).ready(function() {
    hide_plots();
    if (location.pathname.match ("[.][cm]$")) {
	var url = location.pathname.concat("?status");
	$('#status').load(url, function(data) { checkRunning (url, 1000); });
    }
});
