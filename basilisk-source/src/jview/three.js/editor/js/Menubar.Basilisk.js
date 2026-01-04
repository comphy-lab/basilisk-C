import { UIPanel, UIRow, UIText } from './libs/ui.js';

import { BasiliskBufferGeometry } from './BasiliskBufferGeometry.js';

import { AddObjectCommand } from './commands/AddObjectCommand.js';

function MenubarBasilisk( editor ) {
	
	var strings = editor.strings;

	var container = new UIPanel();
	container.setClass( 'menu' );
	container.dom.style.display = 'none';

	var title = new UIPanel();
	title.setClass( 'title' );
	title.setTextContent( strings.getKey( 'menubar/basilisk' ) );
	container.add( title );

	var options = new UIPanel();
	options.setClass( 'options' );
	container.add( options );

	editor.client.signals.oninit.add( function( objects ) {
		for (let type in objects) {
			var option = new UIRow();
			option.setClass( 'option' );
			option.setTextContent( type );
			option.onClick( function () {
				var geometry = new BasiliskBufferGeometry( editor.client, type, objects[type] );
				var mesh = new THREE.Mesh( geometry, new THREE.MeshLambertMaterial( { side: THREE.DoubleSide } ) );
				mesh.name = 'Basilisk';
				BasiliskBufferGeometrySetupObject( mesh, editor );				
				editor.execute( new AddObjectCommand( editor, mesh ) );
			} );
			options.add( option );
		}

		container.dom.style.display = 'block';
	});
	
	editor.client.signals.onclose.add( function() {
		container.dom.style.display = 'none';
		options.clear();
	});
	
	editor.client.signals.onload.add( function( object ) {
		BasiliskBufferGeometrySetupObject( object, editor );
	});
	
	editor.client.signals.onread.add( function( view, geometries ) {
		for ( let geometry of geometries ) {
			var mesh = new THREE.Mesh( geometry, new THREE.MeshLambertMaterial( { side: THREE.DoubleSide } ) );
			mesh.name = 'Basilisk';
			if ( geometry.order )
				mesh.renderOrder = geometry.order;
			BasiliskBufferGeometrySetupObject( mesh, editor );
			editor.scene.add( mesh );
			editor.signals.objectAdded.dispatch( mesh );
		}
		if ( view !== undefined ) {
			editor.camera.position.copy( view.pos );
			editor.camera.quaternion.copy( view.quat );
			editor.signals.cameraResetted.dispatch();
			editor.signals.cameraChanged.dispatch();
		}
		editor.signals.sceneGraphChanged.dispatch();
	});

	return container;
	
}

function BasiliskRead( editor, code )
{
	var read = editor.client.read( code );
	if ( read.view !== undefined ) {
		editor.camera.position.copy( read.view.pos );
		editor.camera.quaternion.copy( read.view.quat );
	}

	for ( let geometry of read.geometries ) {
		var mesh = new THREE.Mesh( geometry, new THREE.MeshLambertMaterial( { side: THREE.DoubleSide } ) );
		mesh.name = 'Basilisk';
		BasiliskBufferGeometrySetupObject( mesh, editor );
		editor.execute( new AddObjectCommand( editor, mesh ) );
	}
	
	editor.signals.cameraResetted.dispatch();
	editor.signals.cameraChanged.dispatch();
}

function sceneHasLight( scene )
{
	var hasLight = false;
	scene.traverse( function ( child ) {		
		if ( child instanceof THREE.Light )
			hasLight = true;	
	});
	return hasLight;
}

function BasiliskBufferGeometrySetupObject( object, editor )
{
	let geometry = object.geometry;

	geometry.parent = object;
	
	geometry.signals.onupdate.add ( function() {

		var updated = geometry.updateParent( geometry.parent );
		if ( updated.parent !== geometry.parent ) {
			var oldObject = geometry.parent;
			var parent = oldObject.parent;
			parent.remove( oldObject );
			geometry.parent = updated.parent;
			parent.add( geometry.parent );
			if ( editor.selected === oldObject )
				editor.select( geometry.parent );
			editor.signals.sceneGraphChanged.dispatch();
		}
		else if ( updated.materialChanged ) {
			editor.signals.materialChanged.dispatch( geometry.parent.material );
			editor.signals.objectChanged.dispatch( geometry.parent );
		}

		if ( geometry.parent.material.type === 'MeshLambertMaterial' && !sceneHasLight( editor.scene ))
			editor.addCameraLight();
		
		editor.signals.geometryUpdated.dispatch( geometry.parent );

		if ( geometry.error === undefined ) {
			if ( geometry.dimension === 2 && editor.camera.dimension === 3 &&
			     editor.client.dimension() === 2 ) {
				editor.camera.dimension = 2;
				editor.camera.position.set( 0, 0, 10 );
				editor.camera.up.set( 0, 1, 0 );
				editor.camera.lookAt( new THREE.Vector3() );
				editor.signals.setDimension.dispatch( 2 );
			}
			else if ( geometry.dimension === 3 && editor.camera.dimension === 2 ) {
				editor.camera.dimension = 3;
				editor.signals.setDimension.dispatch( 3 );
			}
			
			if ( editor.client.refreshed ) {
				editor.client.refreshed = false;
				editor.signals.objectFocused.dispatch( geometry.parent );
			}
		}
		
	});
}

export { MenubarBasilisk, BasiliskRead };
