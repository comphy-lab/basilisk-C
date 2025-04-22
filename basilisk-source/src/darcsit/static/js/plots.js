function hide_plots()
{
    $("[id^=plot]").each(function( index, element ){
	id = $(this).attr('id');
	$('#after' + id).append( $(this) );
	if ($(this).find('#msg_label').length)
	    $('#after' + id).show();
	else
	    $('#after' + id).hide();
	$('#button' + id).click( function(e) {
	    e.preventDefault();
	    $($(this).attr('id').replace('button', '#after')).toggle();
	    return false;
	});
    });
}

/* Call katex to render math */
function render_math()
{
    $('.math').each(function() {
	katex.render($(this).text(), $(this)[0], {
	    throwOnError: false
	});
    });
}
$(function() { render_math(); });
