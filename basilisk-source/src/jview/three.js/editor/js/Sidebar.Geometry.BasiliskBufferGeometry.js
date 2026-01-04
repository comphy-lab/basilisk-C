import { Color } from '../../build/three.module.js';

import { BasiliskBufferGeometry } from './BasiliskBufferGeometry.js';

import { UIRow, UIInput, UIText, UIInteger, UICheckbox, UINumber, UIColor } from './libs/ui.js';

import { SetGeometryCommand } from './commands/SetGeometryCommand.js';

function color_from_rgb(rgb) {
	var decColor = 0x1000000 + Math.floor(rgb[2]*255) + 0x100*Math.floor(rgb[1]*255) + 0x10000*Math.floor(rgb[0]*255);
	return '#' + decColor.toString(16).substr(1);
}

function GeometryParametersPanel( editor, object ) {

	var strings = editor.strings;

	var container = new UIRow();

	var geometry = object.geometry;
	var params = JSON.parse( JSON.stringify( geometry.parameters ) );
	var client = geometry.client;

	var errorRow = new UIRow();
	errorRow.add( new UIText( 'Error' ).setWidth( '90px' ) );
	errorRow.dom.style.color = 'red';
	errorRow.dom.style.display = 'none';
	var error = new UIText( '' ).setWidth( '150px' );
	errorRow.add( error );
	container.add( errorRow );

	if ( geometry.error ) {
		error.setValue( geometry.error );
		errorRow.dom.style.display = 'block';
	}
	else
		errorRow.dom.style.display = 'none';

	var uis = [];

	for (let p in params)
		if (params[p].type) {
		
			var row = new UIRow();
			row.add( new UIText( p ).setWidth( '90px' ) );
			
			var type  = params[p].type;
			var value = params[p].value;
			var ui;
			if (type == 'pstring') {
				ui = new UIInput( value == '(null)' ? '' : value ).setWidth( '150px' );
				ui.setFontSize( '12px' ).onChange(function () {
					params[p].value = this.getValue();
					update();
				} );
			}
			else if (type == 'pbool')
				ui = new UICheckbox( value != false && value != 0 ).onChange( function() {
					params[p].value = this.getValue();
					update();
				} );
			else if (type == 'pint')
				ui = new UIInteger( value ).setWidth( '150px' ).onChange( function () {
					params[p].value = this.getValue();
					update();
				});
			else if (type == 'punsigned')
				ui = new UIInteger( value ).setWidth( '150px' ).onChange( function () {
					params[p].value = this.getValue();
					update();
				}).setRange ( 0, Infinity );
			else if (type == 'pfloat' || type === 'pdouble')
				ui = new UINumber( value ).setPrecision( 8 ).setWidth( '150px' ).onChange( function() {
					params[p].value = this.getValue();
					update();
				});
			else if (type == 'coord') {
				let x = new UINumber( value[0] ).setPrecision( 3 ).setWidth( '50px' ).onChange( function() {
					params[p].value[0] = this.getValue();
					update();
				});
				let y = new UINumber( value[1] ).setPrecision( 3 ).setWidth( '50px' ).onChange( function() {
					params[p].value[1] = this.getValue();
					update();
				});
				ui = new UINumber( value[2] ).setPrecision( 3 ).setWidth( '50px' ).onChange( function() {
					params[p].value[2] = this.getValue();
					update();
				});
				uis.push( x, y );
				row.add( x, y );				
			}
			else if (type == 'color') {
				ui = new UIColor().setValue( color_from_rgb (value) ).onInput( function() {
					const color = new Color( this.getValue() );
					params[p].value[0] = color.r;
					params[p].value[1] = color.g;
					params[p].value[2] = color.b;
					update();
				});
			}

			uis.push( ui );
			row.add( ui );	
			container.add( row );
		}

	if (!editor.client.connected)
		setDisabled( true );	

	function setDisabled( disabled ) {
		for ( let ui of uis )
			ui.dom.disabled = disabled;
	}
	
	function update() {

		object.geometry.dispose();
		var geometry = new BasiliskBufferGeometry( client, params.object, params );
		object.geometry = geometry;
		
		editor.client.signals.onload.dispatch( object );
		editor.signals.historyChanged.dispatch(); // to force saveState without forcing render
	}

	function updated( updatedObject ) {

		if ( updatedObject === object ) {
			if ( object.geometry.error !== undefined ) {
				error.setValue( object.geometry.error );
				errorRow.dom.style.display = 'block';
			}
			else
				errorRow.dom.style.display = 'none';
		}
		
	}
	
	editor.signals.geometryUpdated.add( updated );
	editor.client.signals.onclose.add( function() { setDisabled( true ); } );
	editor.client.signals.oninit.add( function() { setDisabled( false ); } );
	
	return container;

}

export { GeometryParametersPanel };
