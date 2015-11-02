function domain_alias = get_domain()


[~, domain] = system('dnsdomainname');

% is_in_str('lnx.biu.ac.il', domain)
if is_in_str('csail.mit.edu', domain)
    domain_alias = 'csail';
else
    domain_alias = 'cortex';
end

