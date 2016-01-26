// override the default options with something less restrictive.
var options = {
    width: Math.round(parseInt(d3.select('#viewer').style('width'), 10)),
    height: Math.round(parseInt(d3.select('#viewer').style('height'), 10)),
    antialias: true,
    quality: 'high'
};

// insert the viewer under the Dom element with id 'gl'.
var viewer = pv.Viewer(document.getElementById('viewer'), options);

var structure = null;

function loadMolecule() {
    // asynchronously load the PDB file
    // from the server and display it in the viewer.
    $.ajax(window.molecule)
        .done(function(data) {
            structure = pv.io.pdb(data);
            // Match chain colors to D3 representation
            viewer.tube('protein', structure, {
                color: color.byChain(),
                radius: 1.0
            });
            viewer.autoZoom();
        });
}

function setColorForAtom(go, atom, color) {
    var view = go.structure().createEmptyView();
    view.addAtom(atom);
    go.colorBy(pv.color.uniform(color), view);
}

// variable to store the previously picked atom. Required for resetting the color
// whenever the mouse moves.
var prevPicked = null;
// add mouse move event listener to the div element containing the viewer. Whenever
// the mouse moves, use viewer.pick() to get the current atom under the cursor.
parent.addEventListener('mousemove', function(event) {
    var rect = viewer.boundingClientRect();
    var picked = viewer.pick({
        x: event.clientX - rect.left,
        y: event.clientY - rect.top
    });
    if (prevPicked !== null && picked !== null &&
        picked.target() === prevPicked.atom) {
        return;
    }
    if (prevPicked !== null) {
        // reset color of previously picked atom.
        setColorForAtom(prevPicked.node, prevPicked.atom, prevPicked.color);
    }
    if (picked !== null) {
        var atom = picked.target();
        // document.getElementById('picked-atom-name').innerHTML = atom.qualifiedName();
        // get RGBA color and store in the color array, so we know what it was
        // before changing it to the highlight color.
        var color = [0, 0, 0, 0];
        picked.node().getColorForAtom(atom, color);
        prevPicked = {
            atom: atom,
            color: color,
            node: picked.node()
        };

        // highlight chords
        // TODO

        setColorForAtom(picked.node(), atom, 'orange');
    } else {
        // document.getElementById('picked-atom-name').innerHTML = '&nbsp;';
        prevPicked = null;
    }
    viewer.requestRedraw();
});

// load the structure once the DOM has finished loading. That's
// the earliest point the WebGL context is available.
document.addEventListener('DOMContentLoaded', loadMolecule);
