$(function(){
	$('#container>div').not(':first').hide();
	$('#nav>li:first').addClass('active');
	
	var toSection = function(section){
		var link = $('#nav>li>a[href=#'+section+']');
		if(!link){
			return false;
		}
		
		$('#nav>li.active').removeClass('active');
		link.parent().addClass('active');
	
		$('#container>div').hide();
		$('a[name='+section+']').parent().show();
	}
	
	var toHash = function(){
		var hash = location.href.substring(location.href.indexOf("#")+1);
		toSection(hash);
	}
	
	toHash();
	$(window).hashchange(toHash);

	$('textarea')
		.addClass('input-block-level')
		.attr('rows', 10)
		.attr('wrap', 'off')
		.css('font-family', 'Monaco,Menlo,Consolas,"Courier New",monospace');
		
	$('pre.example, #run>pre')
		.addClass('well')
		.addClass('pre-scrollable')
		.css('max-height', '80px');
		
	$('pre.example')
		.each(function(){
			var copy = $('<a href="#">Copy example to input</a>')
								.addClass('pull-right')
								.click(function(){
									$(this).parent().parent().find('textarea').val(
										$(this).parent().next('pre').html().replace(/&gt;/g, '>')
									);
									return false;
								});
								
			$(this).prev().append(copy);
		});
		
	$('#runLink').click(function(){
		$('#run>div.alert').remove();
		var allOK = true;
		$('#run>pre').each(function(){
			var rel = $(this).attr('rel');
			var textarea = $('#'+rel+'>textarea');
			var val = $.trim(textarea.val());
			
			var lines = val.split(/\n/);
			for(var l in lines){
				lines[l] = $.trim(lines[l]).toUpperCase();
			}
			
			lines = lines.filter(function(n){return n}); //clear empty lines
			val = lines.join('\n');
			textarea.val(lines.join('\n'));
			
			var errors = [];
			var toMatch = rel=='msa'
								? /^([-ACDEFGHIKLMNPQRSTVWXY]+)$/
								: /^([-ACDEFGHIKLMNPQRSTVWXY][0-9]+[-ACDEFGHIKLMNPQRSTVWXY])$/;

			for(var l in lines){
				if(!(rel=='msa' && lines[l].substr(0,1)=='>') && !lines[l].match(toMatch)){
					errors.push(parseInt(l)+1);
				}
			}
			
			if((rel=='del' || rel=='neut') && lines.length<2){
				$(this).before('<div class="alert"><strong>WARNING!</strong> At least 2 nsSNPs required</div>');
				allOK = false;
			}
			else if(!lines.length){
				$(this).before('<div class="alert"><strong>WARNING!</strong> Empty</div>');
				allOK = false;
			}
			
			if(errors.length){
				$(this).before('<div class="alert"><strong>WARNING!</strong> Errors were found on lines: '+errors.join(', ') +'</div>');
				allOK = false;
			}
			
			$(this).html(lines.join('\n'));
		});
		
		if(allOK){
			$('#goForIt').removeClass('disabled');
		}
		else {
			$('#goForIt').addClass('disabled');
		}
	});

	$('#goForIt').click(function(){
		if($(this).hasClass('disabled')){
			return false;
		}
		
		$(this)
			.addClass('disabled')
			.html('Running...');

		var toPost = {};
		$('#run>pre').each(function(){
			var rel = $(this).attr('rel');
			var val = $(this).html();
			toPost[rel] = val;
		});

		$.post('./run.php', toPost, function(data){
			$('#goForIt')
				.removeClass('disabled')
				.html('Run');
			$('#results>pre').html(data.response);
			$('#resultLink').tab('show');
		}, 'json');

		return false;
	});
	
});
