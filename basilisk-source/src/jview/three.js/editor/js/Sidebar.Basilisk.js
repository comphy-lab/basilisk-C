import { UIPanel, UIRow, UIText, UIInteger, UIBreak, UIButton, UISpan, UINumber, UIDiv, UIInput } from './libs/ui.js';

class RunButtons {

	constructor( editor ) {
		
		this.runButton = new UIButton( 'Run' ).setWidth( '70px' );
		this.runButton.dom.style['float'] = 'left';
		this.stepButton = new UIButton( 'Step' ).setWidth( '50px' ).setMarginLeft( '7px' );
		this.stepButton.dom.style['float'] = 'left';
		this.state = undefined;
		this.editor = editor;
		var scope = this;		

		this.runButton.onClick( function() { scope.toggleState(); } );

		this.stepButton.onClick( function() {
			scope.editor.client.socket.send( '#Run/Pause:1' );
		});

		this.editor.client.signals.onclose.add( function() { scope.state = undefined; } );
		
		document.addEventListener( 'keydown', function ( event ) {
			if ( scope.state !== undefined && scope.editor.client.connected ) {
				switch ( event.code.toLowerCase() ) {							
				case 'space':
					event.preventDefault();
					scope.toggleState();
					break;
					
				case 'keys': // 's' key
					if ( scope.state === 'run' )
						scope.toggleState();
					else
						scope.editor.client.socket.send( '#Run/Pause:1' );
					break;
				}
			}
		});
	}

	setState( state ) {
		if ( this.state === state )
			return;
		if ( state === 'run' ) {
			this.runButton.dom.textContent = 'Pause';
			this.stepButton.dom.style.display = 'none';
			this.editor.client.socket.send( '#Run/Pause:0' );
			this.editor.client.signals.onrun.dispatch( true );
		}
		else {
			this.runButton.dom.textContent = 'Run';
			this.stepButton.dom.style.display = 'block';
			this.editor.client.socket.send( '#Run/Pause:-1' );
			this.editor.client.signals.onrun.dispatch( false );
		}
		this.state = state;
	}

	toggleState() {
		if ( this.state === undefined )
			return;
		this.setState( this.state === 'run' ? 'pause' : 'run' );
	}
	
}

function SidebarBasilisk( editor ) {

	var strings = editor.strings;
	
	var container = new UIPanel();

	var row = new UIRow();
	row.add( new UIText( 'Basilisk'.toUpperCase() ).setWidth( '90px' ) );
	var host = new UIInput( editor.client.socket.url ).setWidth( '150px' ).setFontSize( '12px' ).onChange( function() {
		if ( editor.client.connected )
			editor.client.socket.close();
		editor.client.socket.url = this.getValue();
		editor.client.socket.open();
	} );
	row.add( host );
	container.add( row );
	
	var connectRow = new UIRow();
	var connect = new UIButton( editor.client.connected ? 'Disconnect' : 'Connect' ).onClick( function () {
		if ( editor.client.connected )
			editor.client.socket.close();
		else
			editor.client.socket.open();
	});
	editor.client.signals.oninit.add( function() {
		connect.setTextContent( 'Disconnect' );
	});
	editor.client.signals.onclose.add( function() {
		connect.setTextContent( 'Connect' );
	});

	connectRow.add( new UIText( 'Server' ).setWidth( '90px' ) );
	connectRow.add( connect );
	container.add( connectRow );

	var div = new UIDiv();
	div.dom.style.display = 'none';
	container.add( div );
	
	var buttons = new RunButtons( editor );
	
	editor.client.signals.oncontrols.add( function( controls ) {
		div.clear();
		for ( let name in controls ) {
			var control = controls[name];
			var row = new UIRow();
			var text = new UIText( name ).setWidth( '90px' );
			row.add( text );
			if ( name === 'Run/Pause' ) {
				buttons.setState( control.value === 0 ? 'run' : 'pause' );
				text.dom.style['float'] = 'left';
				row.add( buttons.runButton, buttons.stepButton );
			}
			else if ( control.type === 'int' ) {
				var control = controls[name];
				var ui = new UIInteger( control.value ).setRange( control.min, control.max ).setWidth( '90px' );
				ui.onChange( function() {
					editor.client.socket.send( '#' + name + ':' + this.getValue() );
				});
				row.add( ui );
			}
			else if ( control.type === 'double' ) {
				var control = controls[name];
				var ui = new UINumber( control.value ).setRange( control.min, control.max ).setWidth( '90px' );
				ui.onChange( function() {
					editor.client.socket.send( '#' + name + ':' + this.getValue() );
				});
				row.add( ui );
			}
			div.add( row );				
		}
		div.dom.style.display = 'block';
	});

	editor.client.signals.onclose.add( function() {
		div.dom.style.display = 'none';
		div.clear();
	});
	
	return container;

}

export { SidebarBasilisk };
