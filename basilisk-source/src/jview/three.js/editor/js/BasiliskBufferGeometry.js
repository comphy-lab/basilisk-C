import { ReconnectingWebSocket } from './libs/reconnecting-websocket.js';

import { BufferGeometry, Float32BufferAttribute, ObjectLoader } from '../../build/three.module.js';

class BasiliskClient {

	constructor( url = 'ws://localhost:7100' ) {

		var Signal = signals.Signal;

		this.socket = new ReconnectingWebSocket( url, null, {
			automaticOpen: false,
			maxReconnectAttempts: 0,
			timeoutInterval: 30000,
			binaryType: 'arraybuffer'
		});
			
		this.geometries = {};
		this.refreshed = true;
		this.objects = [];
		this.controls = {};
		this.signals = {
			oninit:     new Signal(),
			oncontrols: new Signal(),
			onclose:    new Signal(),
			onload:     new Signal(),
			onread:     new Signal(),
			onrun:      new Signal()
		};
		this.connected = false;
		
		var scope = this;

		this.socket.onopen = function() {
			scope.connected = true;
			for (let command in scope.geometries)
				this.send( '+' + command );
		};

		this.socket.onmessage = function( event ) {

			if ( typeof event.data === 'string' ) {
				
				var msg = event.data;
				var first = msg.charAt(0);
				
				if ( first === '{' ) { // the initialisation message from the running Basilisk code

//					console.debug ('[init]', msg);
					scope.objects = JSON.parse (msg);
					scope.signals.oninit.dispatch( scope.objects );

				}
				else if ( first === '#' ) { // the controls from the running Basilisk code

					scope.controls = JSON.parse( msg.substring( 1 ) );
					scope.signals.oncontrols.dispatch( scope.controls );
					
				}
				else { // assume everything else is (text) commands

					if ( msg.substring( 0, 1 ) === '@') { // default display command
						if ( scope.length() > 0 )
							// ignore default display if some geometries are already present
							return;
						msg = msg.substring( 1 );
					}
					var read = scope.read( msg );
					scope.signals.onread.dispatch( read.view, read.geometries );
					
				}
				
			}
			else { /* Assume everything else is binary geometry/app data
			          The format is described in display_send() of [display.h](/src/display.h). */
				
				var bufi = new Uint32Array (event.data, 0, 2);
				var commandlen = bufi[0], errorlen = bufi[1];
				var pos = 8;

				var command = String.fromCharCode.apply(null, new Uint8Array(event.data, pos, commandlen));
				pos += 4*Math.ceil(commandlen/4.);
				
				console.debug ('[command]', command);

				if ( scope.geometries[command] ) {
					var geometries = Object.values (scope.geometries[command]);
					if (errorlen > 0) {
						var error = String.fromCharCode.apply(null, new Uint8Array(event.data, pos, errorlen));
						console.debug ('[error]', error);
						for (let geom of geometries) {
							geom.error = error;
							const indices = [];
							const vertices = [];
							geom.setIndex( indices );
							geom.setAttribute( 'position', new Float32BufferAttribute( vertices, 3 ) );
							geom.deleteAttribute( 'normal' );
							geom.deleteAttribute( 'color' );
							geom.signals.onupdate.dispatch();
						}
					}
					else
						for (let geom of geometries) {
							geom.error = undefined;
							geom.update( event.data, pos );
						}
				}
			}
		};

		this.socket.onclose = function() {
			scope.connected = false;
			scope.signals.onclose.dispatch();
			scope.socket.close();
		};
		
	}
	
	add( geometry ) {

		var command = geometry.command();

		var array = this.geometries[command];
		if ( !array ) {
			this.geometries[command] = array = [];
			if ( this.connected )
				this.socket.send( '+' + command );
		}
			
		if ( !array.includes( geometry ) ) {
			array.push( geometry );
			console.debug ( '[add]', this.geometries );
		}

	}

	remove( geometry ) {

		var command = geometry.command();
		var array = this.geometries[command];
		var index = array !== undefined ? array.indexOf( geometry ) : undefined;
		if ( index != undefined ) {
			array.splice( index, 1 );
			if ( array.length == 0 ) {
				delete this.geometries[command];
				if ( this.connected )
					this.socket.send( '-' + command );
				if ( this.length() === 0 )
					this.refreshed = true;
			}
			console.debug ( '[remove]', this.geometries );
		}
		
	}

	length() {
		
		var n = 0;
		for ( let i in this.geometries )
			n++;
		return n;
		
	}


	dimension() {

		var dimension = 2;
		for ( let command in this.geometries )
			for ( let geometry of this.geometries[command] )
				if ( geometry.dimension > dimension )
					dimension = geometry.dimension;
		return dimension;
		
	}
	
	read( code ) {
		
		function JSONfromCommands( commands ) {
			var s = new String( commands );
			s = s.replace( /([A-Z0-9a-z_]+)[ \t]*=/g, '"$1" :' ).replace( /{/g, '[' ).replace( /}/g, ']' );
			s = s.replace( /'/g, '"' );
			s = s.replace( /([A-Z0-9a-z_]+)[ \t]*\([ \t]*\);/g, '{ "object": "$1" },' );
			s = s.replace( /([A-Z0-9a-z_]+)[ \t]*\(/g, '{ "object": "$1", ' ).replace( /\);/g, ' },' );
			console.debug( s );
			var i = s.lastIndexOf( ',' );
			return JSON.parse( '[' + s.substring(0,i) + ']' + s.substring(i+1) );
		}

		var commands = JSONfromCommands( code );
		var geometries = [];
		var view = undefined;
		for ( let i of commands ) {
			if ( i.object === 'view' ) {
				view = {};
				view.quat = new THREE.Quaternion( i.quat[0], i.quat[1], i.quat[2], i.quat[3] );
				view.pos = new THREE.Vector3( - i.tx, - i.ty, - i.tz );
				view.pos.applyQuaternion( view.quat );
			}
			else {
				var object = this.objects[i.object];
				if ( object ) { // the Basilisk client knows this object type
					var command = i.object + ' (', sep = '';
					var type = i.object;
					delete i.object;
					for ( let j in i ) {
						command += sep + j + ' = ';
						if ( typeof(i[j]) === 'string')
							command += '"' + i[j] + '"';
						else
							command += i[j];
						sep = ', ';
					}
					command += ');';
					if ( this.geometries[command] === undefined ) {
						// this is a new Basilisk geometry
						var params = JSON.parse( JSON.stringify( object ) );
						var error = false, order = undefined;
						for ( let j in i ) {
							if ( j === 'order' )
								order = i[j];
							else if ( params[j] === undefined ) {
								alert( "Invalid parameter '" + j + "' in '" + command + "'" );
								error = true;
							}
							else
								params[j].value = i[j];
						}
						if ( !error ) {
							var geometry = new BasiliskBufferGeometry( this, type, params );
							if ( order !== undefined )
								geometry.order = order;
							geometries.push( geometry );
						}
					}
				}
			}
		}

		return { geometries: geometries, view: view };
	}
}

class BasiliskBufferGeometry extends BufferGeometry {

	constructor( client, type, parameters ) {

		super();
		this.type = 'BasiliskBufferGeometry';

		this.client = client;

		var Signal = signals.Signal;

		this.signals = {
			onupdate:   new Signal()
		};
		
		this.parameters = JSON.parse( JSON.stringify( parameters ) );
		this.parameters.object = type;
		this.name = type;

		const indices = [];
		const vertices = [];
		
		this.setIndex( indices );
		this.setAttribute( 'position', new Float32BufferAttribute( vertices, 3 ) );		
		this.client.add( this );

	}

	command() {
		
		var s = this.parameters.object + ' (';
		var sep = '';
		for (let p in this.parameters) {
			let type  = this.parameters[p].type;
			let value = this.parameters[p].value;
			if (type == 'pstring') {
				if (value != '(null)' && value != '') {
					s += sep + p + ' = "' + value + '"';
					sep = ', ';
				}
			}
			else if (type == 'pbool') {
				if (value != '0') {
					s += sep + p + ' = true';
					sep = ', ';
				}			
			}
			else if (type == 'pint' || type == 'punsigned') {
				value = parseInt (value);
				if (value != 0) {
					s += sep + p + ' = ' + value;
					sep = ', ';
				}
			}
			else if (type == 'pfloat' || type === 'pdouble') {
				value = parseFloat (value);
				if (value != 0) {
					s += sep + p + ' = ' + value;
					sep = ', ';
				}			
			}
			else if (type == 'coord' || type == 'color') {
				if (value[0] != 0 || value[1] != 0 || value[2] != 0) {
					s += sep + p + ' = {' + value[0] + ',' + value[1] + ',' + value[2] + '}';
					sep = ', ';				
				}
			}
		}
		s += ');'
		return s;

	}

	update( buffer, pos ) {
			
		var bufi = new Uint32Array( buffer, pos, 6 );
		pos += 24;
		var dimension = bufi[0];
		var type = bufi[1], positionlen = bufi[2], normallen = bufi[3], colorlen = bufi[4], indexlen = bufi[5];

		this.dimension = dimension;
		this.parentType = type;		

		console.debug ('[geometry]', positionlen, normallen, colorlen, indexlen, type);
		
		var positions = new Float32Array (buffer, pos, positionlen/4);
		pos += positionlen;
		this.setAttribute( 'position', new THREE.BufferAttribute( positions, 3, false ) );

		if ( normallen > 0 ) {
			var normals = new Float32Array (buffer, pos, normallen/4);
			pos += normallen;
			this.setAttribute( 'normal', new THREE.BufferAttribute( normals, 3, false ) );
		}
		else
			this.deleteAttribute( 'normal' );

		if ( colorlen > 0 ) {
			var colors = colorlen >= positionlen ?
			    new Float32Array (buffer, pos + colorlen - positionlen, positionlen/4) :
			    new Float32Array (buffer, pos, colorlen/4);
			pos += colorlen;
			this.setAttribute( 'color', new THREE.BufferAttribute( colors, 3, false ) );
		}
		else
			this.deleteAttribute( 'color' );

		if ( indexlen > 0 ) {
			var indices = new Uint32Array (buffer, pos, indexlen/4);
			pos += indexlen;
			this.setIndex( new THREE.BufferAttribute( indices, 1 ) );
		}
		else {
			var indices = [];
			this.setIndex( indices );
		}
		
		this.computeBoundingBox();
		this.computeBoundingSphere();

		this.signals.onupdate.dispatch();
	}

	setParentColor( parent ) {
		
		if ( this.attributes.color && this.attributes.color.count === 1 ) {
			var colors =  this.attributes.color.array;
			var c = new THREE.Color (colors[0], colors[1], colors[2]);
			if ( !parent.material.color.equals(c) ) {
				parent.material.color.copy (c);
				return true;
			}
		}
		
		return false;
		
	}
	
	updateParent( parent ) {

		var hasColor = this.attributes.color && this.attributes.color.count === this.attributes.position.count ? true : false;
		
		if ( this.parentType === 0 && parent.type !== 'LineSegments' ) {
			const material = new THREE.LineBasicMaterial( { vertexColors: hasColor } );
			var newParent = new THREE.LineSegments( this, material );
			newParent.name = parent.name;
			newParent.renderOrder = 1;
			this.setParentColor( newParent );
			return { parent: newParent, materialChanged: false };
		}
		else if ( this.parentType === 1 && parent.type !== 'Mesh' ) {
			var newParent =
			    new THREE.Mesh (this, this.hasAttribute( 'normal' ) ?
					    new THREE.MeshLambertMaterial( { side: THREE.DoubleSide, vertexColors: hasColor }) :
					    new THREE.MeshBasicMaterial ({ side: THREE.DoubleSide, vertexColors: hasColor }));
			newParent.name = parent.name;
			newParent.renderOrder = 0;
			this.setParentColor( newParent );
			return { parent: newParent, materialChanged: false };
		}
		
		var changed = false;
		var hasNormal = this.attributes.normal && this.attributes.normal.count === this.attributes.position.count ? true : false;
		if ( this.parentType == 1 && !hasNormal && parent.material.type !== 'MeshBasicMaterial' ) {
			parent.material.dispose();
			parent.material = new THREE.MeshBasicMaterial ({ side: THREE.DoubleSide, vertexColors: hasColor });
			changed = true;
		}
		
		if ( hasColor ) {
			if ( parent.material.vertexColors == false && parent.material.color != undefined ) {
				parent.material.vertexColors = true;
				parent.material.color.setRGB (1, 1, 1);
				parent.material.needsUpdate = true;
				changed = true;
			}
		}
		else if ( parent.material.vertexColors == true ) {
			parent.material.vertexColors = false;
			parent.material.needsUpdate = true;
			changed = true;
		}

		if ( this.setParentColor( parent ) )
			changed = true;
		
		return { parent: parent, materialChanged: changed };

	}

	materialName() {

		for ( let p in this.parameters ) {
			let param = this.parameters[p];
			if ( param.type === 'pstring' && param.value !== '' && param.value !== '(null)' )
				return param.value;
		}
		return '';
		
	}
	
	dispose() {

		this.client.remove( this );
		super.dispose();
		
	}

}

function BasiliskGeometryFromJSON( json, client )
{
	var geometry = new BasiliskBufferGeometry( client, json.object, json );
	var parameters = geometry.parameters;
	for ( let p in parameters )
		if ( p !== 'object' && ( parameters[p].type === undefined ||
					 parameters[p].cardinality === undefined ) )
			delete parameters[p];
	
	geometry.uuid = json.uuid;
	if ( json.name !== undefined ) geometry.name = json.name;
	if ( geometry.isBufferGeometry === true && json.userData !== undefined )
		geometry.userData = data.userData;

	return geometry;
}

class BasiliskObjectLoader extends ObjectLoader {

	constructor( client, manager ) {

		super( manager );		
		this.client = client;
		
	}

	parseGeometries( json, shapes ) {

		const geometries = {};
		
		if ( json !== undefined ) {
		
			for ( let i = 0; i < json.length; i ++ ) {

				const data = json[ i ];

				if ( data.type === 'BasiliskBufferGeometry' )  {

					geometries[ data.uuid ] = BasiliskGeometryFromJSON( data, this.client );
					json.splice( i, 1 ); i--;
				
				}
				
			}
			
		}

		return Object.assign( super.parseGeometries( json, shapes ), geometries );

	}
	
	parseObject( data, geometries, materials, animations ) {
		
		let object = super.parseObject( data, geometries, materials, animations );
		if ( object.geometry && object.geometry.type === 'BasiliskBufferGeometry' )
			this.client.signals.onload.dispatch( object );
		
		return object;
		
	}

}

export { BasiliskClient, BasiliskBufferGeometry, BasiliskObjectLoader, BasiliskGeometryFromJSON };
