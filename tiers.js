var EXAMPLES = {};

var MODIFY = {
  "dolgo.1" : "DOLGO / S_",
  "dolgo.2" : "DOLGO / _S",
  "sca.1" : "SCA / S_",
  "sca.2" : "SCA / _S",
  "cv.1" : "SCA / S_",
  "cv.2" : "SCA / _S",
  "asjp.1" : "ASJP / S_",
  "asjp.2" : "ASJP / _S",
  "word_length.0" : "LENGTH",
  "position.0" : "POSITION",
  "prostring.0" : "PROSODY"
};

function savePlotWord(sound, tag, classes) {

  if (sound in DOLGO || sound == '-') {
    return plotWord(sound, tag, classes);
  }
  else if (sound[0] in DOLGO) {
    return plotWord(sound, tag, classes);
  }
  else {
    return sound;
  }
}

function show_data (sound, idx) {
  /* present the data in tabular form */

  if (typeof idx == 'undefined') {
    idx = -1;
  }
  
  var visited = [];
  var examples = {};

  var basic_taxon = '';
  var row_color = 'white';
  var basic_reflex = '';


  var subset = [];
  var active_tiers = [];
  for (var i=0, row; row=PATTERNS[i]; i++) {
    if (row[1] == sound || sound=='all') {
      subset.push(row);
      for (tier in row[3]) {
	active_tiers.push(tier);
      }
    }
  }
  var columnDefs = [
    {targets:[0],orderData: [0]},
    {targets:[1],orderData: [0,1]},
    {targets:[2],orderData: [0,2]}
  ];

  var txt = '<table style="overflow-x:hidden;width:60%;margin:5px;" id="datatable" border="1"><thead>';
  txt += '<tr><th>LANGUAGE</th><th >SOURCE</th><th>TARGET</th>';
  var number_of_tiers = 0;
  var mycounter = 3
  for (var i=0,tier; tier=TIERS[i]; i++) {
    if (active_tiers.indexOf(tier) != -1) {
      number_of_tiers += 1;
      txt += '<th>'+MODIFY[tier]+'</th>';
      columnDefs.push({targets:[mycounter],orderData:[0,1,2,mycounter]});
      mycounter += 1;
    }
  }
  txt += '<th>EXAMPLES</th></tr></thead><trows>';

  var table = [];

  for (var i=0,row; row=subset[i]; i++) {
    var taxon = row[0];
    var source = row[1];
    var target = row[2];
    var context = row[3];
    var proto = row[4];
    var reflex = row[5];
    var pos = row[6];
    var cogid = row[7];

    if (taxon != basic_taxon) {
      basic_taxon = taxon;
      //txt += '<tr><td style="background:black;border:1px solid black;" colspan="'+(number_of_tiers+4)+'"></td></tr>';
      row_color = "white";

    }
    if (basic_reflex != target) {
      basic_reflex = target;
      //var row_style="border: 2px solid black;";
      //txt += '<tr><td style="border:1px solid black;" colspan="'+(number_of_tiers+4)+'"></td></tr>';
      if (row_color == 'white') {
	row_color='white';
      }
      else {
	row_color = 'white';
      }
    }
    else {
      var row_style="";
    }


    var check = taxon+'/'+source+'/'+target;
    for (key in context) {
      check += '/'+key+':'+context[key];
    }
    
    if (visited.indexOf(check) == -1) {
      visited.push(check);
      txt += '<tr style="background-color:'+row_color+'">';
    	txt += '<td>'+taxon+'</td><td>'+savePlotWord(source, 'span')+'</td><td>'+savePlotWord(target, 'span')+'</td>';
      for (var j=0,tier; tier=TIERS[j]; j++) {
        if (active_tiers.indexOf(tier) != -1) {
          if (tier in context) {
            txt += '<td>'+savePlotWord(context[tier], 'span')+'</td>';

          }
          else {
            txt += '<td></td>';
          }
        }
      }
      txt += '<td id="'+check+'" style="cursor:pointer" onclick="showAlignments(\''+check+'\');"></td>';
      examples[check] = [[proto, reflex, pos, cogid]];

      txt += '</tr>';
    }
    else {
      examples[check].push([proto, reflex, pos, cogid])
    }
  }

  txt += '</trows></table>';
  columnDefs.push({targets: [mycounter], orderData:[0,1,2,mycounter]});

  document.getElementById('data').innerHTML = txt;
  for (key in examples) {
    document.getElementById(key).innerHTML = examples[key].length;
  }
  
  $('#datatable').DataTable({scrollCollapse:true,columnDefs:columnDefs,scrollY: "50vh", paging:false});

  EXAMPLES = examples;
  console.log(examples);

}

function showAlignments(key) {
  
  var txt = '<div style="justify-content: space-between;display:flex; flex-wrap:wrap;" id="alignment_table">';
  var counter = 1;
  for (var i=0,row; row=EXAMPLES[key][i]; i++) {
    var alm1 = row[0];
    var alm2 = row[1];
    pos = row[2];
    cogid = row[3];
  
    txt += '<table border="1" style="margin:5px;"><tr><th>Cognate Set '+cogid+'</th></tr><tr><td>';
    counter += 1;
    for (var j=0,sound; sound=alm1[j]; j++) {
      if (j == pos) {
	txt += savePlotWord(sound, 'span', 'alert');
      }
      else {
	txt += savePlotWord(sound, 'span');
      }
    }
    txt += '</td></tr>';
    txt += '<tr><td>';
    for (var j=0,sound; sound=alm2[j]; j++) {
      
      if (j == pos) {
	txt += savePlotWord(sound, 'span', 'alert');
      }
      else {
	txt += savePlotWord(sound, 'span');
      }
    }
    txt += '</td></tr></table>';
  }
  txt += '</div>';
  document.getElementById('alignments').innerHTML = txt;
}
function prepare_data() {
  var txt = '<table id="sound_table"><tr>';
  count = 0;
  for (var i=0,sound; sound=SOUNDS[i]; i++) {

    if (count % 8 == 0) {
      txt += '</tr>';
      if (i+1 != SOUNDS.length) {
	txt += '<tr>';
      }
    }

    txt += '<td onclick="show_data(\''+sound+'\')">'+savePlotWord(sound)+'</td>';
    count += 1;
  }
  txt += '</tr></table>';
  document.getElementById('sounds').innerHTML = txt;
}

prepare_data()
show_data('p');

