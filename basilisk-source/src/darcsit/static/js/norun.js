function run() {
    $('#status').html("<div class=error><div id=msg_logo><img src=/img/error.png></div><div id=msg_label>Please ask the system administrator for permission to run code on the server.</div></div>");
};

$(document).ready(function(){
    $("#runButton").show();
});
