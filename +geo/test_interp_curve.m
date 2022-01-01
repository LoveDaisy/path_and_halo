function test_interp_curve()
test_cases = {@suite1};
num = length(test_cases);

fprintf('Start testing for interp_curve...\n');
for i = 1:num
    fprintf('Testing suite %d/%d...\n', i, num);
    test_cases{i}();
    fprintf('suite %d/%d passed!\n', i, num);
end
end

function suite1()
pts0 = [0, 0;
    1, 0;
    1, 1;
    0, 1];
ds = 0.1;

pts_e = [0, 0;
    0.121204671510824, -0.109902568131678;
    0.235290625335865, -0.195414615624755;
    0.342257886832422, -0.258011959258022;
    0.442106481357790, -0.299170415810269;
    0.534836434269266, -0.320365802060285;
    0.620447770924148, -0.323073934786863;
    0.698940516679732, -0.308770630768792;
    0.770314696893315, -0.278931706784862;
    0.834570336922195, -0.235032979613864;
    0.891707462123667, -0.178550266034588;
    0.941726097855029, -0.110959382825825;
    0.984626269473578, -0.0337361467663645;
    1, 0;
    1.02040800233661, 0.0516436253650022;
    1.04907132180142, 0.143704116789485;
    1.07061625322531, 0.240969510728294;
    1.08504282196558, 0.341963990402639;
    1.09235105337951, 0.445211739033728;
    1.09254097282442, 0.549236939842772;
    1.08561260565759, 0.652563776050980;
    1.07156597723632, 0.753716430879562;
    1.05040111291791, 0.851219087549728;
    1.02211803805965, 0.943595929282687;
    1, 1;
    0.986716778018851, 1.02937113929965;
    0.944197358152798, 1.10706890082182;
    0.894559803818791, 1.17521339707042;
    0.837804140374128, 1.23232881126665;
    0.773930393176105, 1.27693932663172;
    0.702938587582018, 1.30756912638684;
    0.624828748949166, 1.32274239375322;
    0.539600902634844, 1.32098331195207;
    0.447255073996350, 1.30081606420460;
    0.347791288390980, 1.26076483373201;
    0.241209571176032, 1.19935380375553;
    0.127509947708802, 1.11510715749636;
    0.00669244334658725, 1.00654907817570;
    0, 1];

s_e = [0;
    0.163612795581473;
    0.306188793304188;
    0.430125964733965;
    0.538124855350742;
    0.633246296494430;
    0.718900455451834;
    0.798685765140298;
    0.876046187745658;
    0.953865758052584;
    1.03420844913004;
    1.11829406866969;
    1.20663348746106;
    1.24370745661653;
    1.29923718379226;
    1.39565668370830;
    1.49527967768124;
    1.59729933662161;
    1.70080541298328;
    1.80483078716073;
    1.90838964652291;
    2.01051294136241;
    2.11028628218729;
    2.20689587390567;
    2.26748157189697;
    2.29971676710282;
    2.38828788661309;
    2.47259422782914;
    2.55311364084147;
    2.63102354818485;
    2.70834121172008;
    2.78791114909483;
    2.87315714698424;
    2.96767947211771;
    3.07490421604056;
    3.19791225255549;
    3.33942232916669;
    3.50184685301116;
    3.51121057144450];

[pts, s] = geo.interp_curve(pts0, ds);

assert(all(abs(pts(:) - pts_e(:)) < 1e-6));
assert(all(abs(s - s_e) < 1e-6));
end
