#!/usr/bin/env node

var fs = require("fs");
const katex = require('katex');
const jsdom = require("jsdom");
const { JSDOM } = jsdom;

const argv = require('yargs')
  .scriptName("katex_serial")
  .usage("$0 -i inputfile -o outputfile")
  .option('input',{
    alias: 'i',
    describe: 'input file to process KaTeX maths'
  })
  .option('output',{
    alias: 'o',
    describe: 'output file with KaTeX rendered'
  })
  .demandOption(['input','output'],"Please provide both input and output arguments to work this tool")
  .help()
  .alias('help','h')
  .argv

JSDOM.fromFile( argv['input'] ,{}).then(dom => {
  var mathElements = dom.window.document.getElementsByClassName("math");
  var macros = [];
  for (var i = 0; i < mathElements.length; i++) {
    var texText = mathElements[i].firstChild;
    if (mathElements[i].tagName == "SPAN") {
      mathElements[i].innerHTML = katex.renderToString(texText.data, {
        displayMode: mathElements[i].classList.contains('display'),
        throwOnError: false,
        macros: macros,
        fleqn: false
      });
    }
    
  }

  fs.writeFile(argv['output'], dom.serialize(), (err) => {
    if (err) console.log(err);
    console.log("Successfully KaTeX rended!");
  });

});
