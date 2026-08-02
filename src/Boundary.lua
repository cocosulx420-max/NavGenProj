--!strict
-- NavGen.Boundary — clean outlines from the surfel grid alone.
--
-- WHAT PROBLEM THIS SOLVES, PRECISELY. The floor stage raycasts DOWN onto the
-- real part for every cell, so heights are already exact: ramps come out smooth
-- and stair risers stay crisp. There is no vertical staircase in this pipeline.
-- What is left is the XZ OUTLINE. Cells are axis-aligned and walls are not, so
-- a rotated wall leaves a jagged run of cells along its base. That jaggedness is
-- what inflates polygon counts and manufactures junk portals.
--
-- THE RULE THIS MODULE OBEYS: it never asks a Part anything. No raycasts, no
-- overlap queries, no CFrames, no face planes. The only contact with real
-- geometry in the whole chain is the raycast the floor stage already paid for.
-- Everything here is arithmetic on the surfel grid.
--
-- That is the entire reason it survives real maps. The previous approach needed
-- a readable planar face with a usable CFrame to steal a line from — which a
-- Union or a MeshPart does not have, and which interpenetrating parts corrupt.
-- A raycast does not care what it hit.
--
-- THE STAIRCASE DIES IN STEP 4. A best-fit line is grown along the boundary
-- cells while the MAXIMUM perpendicular residual stays under a stud. The cells
-- of a rotated wall's staircase all sit within a stud of one straight line at
-- the true angle, so they collapse into a single segment and the fit recovers
-- the angle without ever being told it. Corners are then a BYPRODUCT — they are
-- where the fit fails — rather than something detected separately against real
-- geometry, which is what every previous attempt died of.
--
-- SAFETY PROPERTY, and it is structural rather than tuned: the only morphology
-- here is EROSION. Cells are never dilated, and the offset only ever moves a
-- line inward. Erosion can remove walkable ground but it can never invent
-- connectivity, so no amount of it can weld two rooms through a wall. A
-- dilate-then-erode close would bridge any wall thinner than its kernel, which
-- is the through-wall portal bug manufactured on purpose. It is not used here
-- and must not be added.

local Boundary = {}

export type Config = {
	-- Layer separation. Two adjacent cells are the same layer when their height
	-- difference is under this. It must EXCEED the rise across one cell on the
	-- steepest walkable slope -- tan(65 deg) = 2.14 -- or every steep ramp
	-- shatters into one layer per cell.
	stepTol: number?,
	minClearance: number?,

	-- Greedy fit. The maximum perpendicular distance a boundary cell centre may
	-- sit from its segment's line. One cell: anything the raster can express as
	-- "straight" is straight.
	fitTol: number?,

	-- Inward offset. Pushed by clamp(maxD - margin, 0, agentRadius), GRADED
	-- rather than switched: a hard "skip the offset when the corridor is narrow"
	-- makes two adjacent polygons straddling the threshold offset by r and by 0,
	-- so they no longer share an edge. Grading is continuous, and narrow
	-- corridors keep their walkable area instead of vanishing.
	agentRadius: number?,
	margin: number?,

	-- Cleanup. Merging happens BEFORE corners are intersected: two very short
	-- adjacent segments are nearly parallel, and their intersection is then
	-- numerically unstable.
	minSegLen: number?,
	collinearDeg: number?,

	-- Acute corners send two offset lines' intersection arbitrarily far out --
	-- the classic miter spike. Past this distance from the corner it replaces,
	-- bevel instead.
	miterLimit: number?,
}

local DEFAULT = {
	stepTol = 2.2,
	minClearance = 1.5,
	fitTol = 1.0,
	agentRadius = 1.5,
	margin = 0.5,
	minSegLen = 1.0,
	collinearDeg = 5,
	miterLimit = 3.0,
}

local function merged(cfg): any
	local c = {}
	for k, v in pairs(DEFAULT) do c[k] = v end
	if cfg then for k, v in pairs(cfg) do if v ~= nil then c[k] = v end end end
	return c
end

local INF = 1e20

--------------------------------------------------------------------------
-- STEP 2 — exact Euclidean distance transform (Felzenszwalb & Huttenlocher)
--
-- It must be TRUE Euclidean. Stepping over 4-neighbours gives a diamond kernel
-- and 8-neighbours a square one, so an offset derived from either comes out
-- wrong on the diagonals -- which is exactly where the rotated walls are.
--------------------------------------------------------------------------

local function edt1d(f: {number}, n: number): {number}
	local v, z, d = table.create(n, 0), table.create(n + 1, 0), table.create(n, 0)
	local k = 1
	v[1] = 1; z[1] = -INF; z[2] = INF
	for q = 2, n do
		local s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k])
		while s <= z[k] do
			k -= 1
			s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k])
		end
		k += 1; v[k] = q; z[k] = s; z[k + 1] = INF
	end
	k = 1
	for q = 1, n do
		while z[k + 1] < q do k += 1 end
		local dq = q - v[k]
		d[q] = dq * dq + f[v[k]]
	end
	return d
end

-- `inside(i, j)` over a w x h raster -> squared distance to the nearest cell
-- that is NOT inside.
local function edt2d(w: number, h: number, inside: (number, number) -> boolean): {number}
	local g = table.create(w * h, 0)
	local col = table.create(h, 0)
	for i = 1, w do
		for j = 1, h do col[j] = inside(i, j) and INF or 0 end
		local d = edt1d(col, h)
		for j = 1, h do g[(j - 1) * w + i] = d[j] end
	end
	local row = table.create(w, 0)
	for j = 1, h do
		local base = (j - 1) * w
		for i = 1, w do row[i] = g[base + i] end
		local d = edt1d(row, w)
		for i = 1, w do g[base + i] = d[i] end
	end
	return g
end

--------------------------------------------------------------------------
-- small 2D helpers, all in world XZ
--------------------------------------------------------------------------

type P2 = { x: number, z: number }

local function sub(a: P2, b: P2): P2 return { x = a.x - b.x, z = a.z - b.z } end
local function dot(a: P2, b: P2): number return a.x * b.x + a.z * b.z end
local function len(a: P2): number return math.sqrt(a.x * a.x + a.z * a.z) end

-- Interior lies to the LEFT of travel, so this is the outward side.
local function outwardOf(d: P2): P2
	local m = len(d)
	if m < 1e-9 then return { x = 0, z = 0 } end
	return { x = d.z / m, z = -d.x / m }
end

-- Total least squares: the principal axis of the point set. NOT ordinary least
-- squares -- a boundary run can be near-vertical in XZ, where fitting z as a
-- function of x blows up.
--
-- It takes an explicit list of point INDICES rather than a range: the run that
-- straddles the loop's arbitrary start point wraps, and a range cannot express
-- that.
local function fitLine(pts: {P2}, idx: {number}): (P2, P2)
	local n = #idx
	local sx, sz = 0, 0
	for _, i in ipairs(idx) do sx += pts[i].x; sz += pts[i].z end
	local cx, cz = sx / n, sz / n
	local sxx, szz, sxz = 0, 0, 0
	for _, i in ipairs(idx) do
		local dx, dz = pts[i].x - cx, pts[i].z - cz
		sxx += dx * dx; szz += dz * dz; sxz += dx * dz
	end
	-- larger eigenvalue of the 2x2 covariance, closed form
	local tr, det = sxx + szz, sxx * szz - sxz * sxz
	local disc = math.max(tr * tr * 0.25 - det, 0)
	local lam = tr * 0.5 + math.sqrt(disc)
	local dx, dz
	if math.abs(sxz) > 1e-12 then
		dx, dz = lam - szz, sxz
	elseif sxx >= szz then
		dx, dz = 1, 0
	else
		dx, dz = 0, 1
	end
	-- ORIENT ALONG TRAVEL. A principal axis has no sign -- PCA is just as happy
	-- to hand back the reverse direction -- but everything downstream reads the
	-- sign as meaning something. outwardOf() assumes interior-on-the-left, which
	-- is only true if dir runs the way the loop was walked; a flipped segment
	-- gets its "outward" normal pointing INTO the floor, so step 5 biases the
	-- wrong way and step 6 offsets outward, straight through the wall it was
	-- supposed to stand off from. Point it along the run's own span.
	local sxs = pts[idx[n]].x - pts[idx[1]].x
	local szs = pts[idx[n]].z - pts[idx[1]].z
	if dx * sxs + dz * szs < 0 then dx, dz = -dx, -dz end
	local m = math.sqrt(dx * dx + dz * dz)
	if m < 1e-12 then dx, dz, m = 1, 0, 1 end
	return { x = cx, z = cz }, { x = dx / m, z = dz / m }
end

local function maxResidual(pts: {P2}, idx: {number}, c: P2, d: P2): number
	local nx, nz = -d.z, d.x
	local worst = 0
	for _, i in ipairs(idx) do
		local r = math.abs((pts[i].x - c.x) * nx + (pts[i].z - c.z) * nz)
		if r > worst then worst = r end
	end
	return worst
end

--------------------------------------------------------------------------
-- STEP 1 — layers by connectivity AND 2D injectivity, never by height band
--
-- A height slice would put a balcony and the courtyard beneath it in the same
-- layer and corrupt everything downstream, so connectivity is the right primary
-- relation. But connectivity ALONE is not enough, and the reason is not exotic.
--
-- The claim it replaced was that a balcony over a courtyard is safe because the
-- two are separate components. That is true only while they are disconnected.
-- Put a ramp between them -- the ordinary case, and the case this pipeline
-- exists to handle -- and they are one component, which then gets flattened
-- into ONE 2D raster where both floors compete for the same "x:z" key. The
-- loser is silently dropped and every stage downstream reads a floor plan that
-- is part ground and part balcony.
--
-- That is not a corner case. On SmallMap it was 17.9% of the largest layer,
-- every one of those cells holding heights more than 5 studs apart and 6243 of
-- them more than 10 apart, worst pair 36.5 studs.
--
-- So a layer is grown under BOTH constraints: a cell may join only if it is
-- adjacent within the step tolerance AND the layer does not already occupy that
-- x:z at an incompatible height. Cells refused on the second test are not lost;
-- they seed the next layer, which is exactly the balcony peeling off the
-- courtyard. Every layer is now injective on x:z by construction, which is the
-- precondition the 2D raster in steps 2-8 always silently assumed.
--------------------------------------------------------------------------

function Boundary.layers(floorData: any, cfg: Config?)
	local c = merged(cfg)

	local nodes: {any} = {}
	local at: {[string]: {number}} = {}
	-- Every surface height present at each cell, INCLUDING ones filtered out of
	-- the walkable set. It is what lets a boundary edge be classified without
	-- asking a Part anything: a wall is a column whose surface is above us, a
	-- dropoff is one whose surface is below.
	local allY: {[string]: {number}} = {}
	for k, bucket in pairs(floorData.index) do
		local ys = {}
		for _, s in ipairs(bucket) do ys[#ys + 1] = s.pos.Y end
		allY[k] = ys

		-- A ClipRamp is the smooth surface authored over a staircase, so where
		-- one covers a cell it IS the floor and the risers underneath are not.
		-- Keeping both makes every step its own micro-layer, which is what drew
		-- the picket fences over both staircases.
		local hasClip = false
		for _, s in ipairs(bucket) do
			if s.clip then hasClip = true; break end
		end

		for _, s in ipairs(bucket) do
			if s.clearance >= c.minClearance and not (hasClip and not s.clip) then
				local p = s.pos
				local ix, iz = math.floor(p.X), math.floor(p.Z)
				nodes[#nodes + 1] = { ix = ix, iz = iz, y = p.Y, surfel = s, comp = 0, clip = s.clip }
				local b = at[k]
				if not b then b = {}; at[k] = b end
				b[#b + 1] = #nodes
			end
		end
	end

	local DIRS = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } }
	local comps: {any} = {}
	for i = 1, #nodes do
		if nodes[i].comp == 0 then
			local id = #comps + 1
			local members = { i }
			nodes[i].comp = id
			-- the height this layer has claimed at each x:z. Its existence is what
			-- keeps the layer injective, and injectivity is what makes the 2D
			-- raster downstream mean anything.
			local claimed: {[string]: number} = {}
			claimed[nodes[i].ix .. ":" .. nodes[i].iz] = nodes[i].y
			local queue, head = { i }, 1
			while head <= #queue do
				local cur = nodes[queue[head]]; head += 1
				for _, d in ipairs(DIRS) do
					local key = (cur.ix + d[1]) .. ":" .. (cur.iz + d[2])
					local b = at[key]
					if b then
						for _, ni in ipairs(b) do
							local nb = nodes[ni]
							if nb.comp == 0 and math.abs(nb.y - cur.y) <= c.stepTol then
								-- Refuse a cell this layer already holds at another
								-- height. Refusing is not discarding: nb keeps comp 0
								-- and seeds a later layer, which is the balcony
								-- separating itself from the courtyard it is ramped to.
								local held = claimed[key]
								if not held or math.abs(held - nb.y) <= c.stepTol then
									nb.comp = id
									if not held then claimed[key] = nb.y end
									members[#members + 1] = ni
									queue[#queue + 1] = ni
								end
							end
						end
					end
				end
			end
			comps[id] = members
		end
	end

	return nodes, comps, allY
end

--------------------------------------------------------------------------
-- STEP 3 — trace the contour
--
-- From where SURFELS end, never from SVO solid voxels: the octree is
-- deliberately conservative and would inflate the outline by up to a leaf.
--
-- Walls and cliff edges are captured identically -- both are just "floor stops
-- here" -- so there is no special case for either.
--
-- Boundary EDGES are chained rather than boundary cells walked, because the
-- boundary of a set of cells is a closed loop by construction: there is no open
-- end to chase and no tolerance involved. The ordered cells for the fit fall
-- out of the edge order.
--------------------------------------------------------------------------

local function traceLoops(member: {[string]: boolean}, cells: {any}, cellY: {[string]: number}, stepTol: number)
	local DIRS = { { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } }
	-- lattice corners of cell (x,z)'s edge in direction di, oriented so the loop
	-- runs with the interior on its left
	local function corners(x: number, z: number, di: number)
		if di == 1 then return { x + 1, z }, { x + 1, z + 1 }
		elseif di == 2 then return { x + 1, z + 1 }, { x, z + 1 }
		elseif di == 3 then return { x, z + 1 }, { x, z }
		else return { x, z }, { x + 1, z } end
	end

	local segs: {any} = {}
	local byStart: {[string]: {any}} = {}
	for _, cell in ipairs(cells) do
		for di, d in ipairs(DIRS) do
			local nk = (cell.ix + d[1]) .. ":" .. (cell.iz + d[2])
			-- Floor stops here if there is no neighbour OR if the neighbour is a
			-- cliff away. Membership alone is not the test: a layer is connected
			-- in 3D, but its 2D projection is not, so a ramp climbing beside the
			-- floor it started from puts both in one layer with a 15-stud drop
			-- between cells 1.7 studs apart. Without this the raster reads that
			-- drop as continuous ground and the polygon spans the cliff.
			local drop = cellY[nk]
			if not member[nk] or (drop ~= nil and math.abs(drop - cell.y) > stepTol) then
				local a, b = corners(cell.ix, cell.iz, di)
				local s = { a = a, b = b, cell = cell, nk = nk }
				segs[#segs + 1] = s
				local k = a[1] .. "," .. a[2]
				local bucket = byStart[k]
				if not bucket then bucket = {}; byStart[k] = bucket end
				bucket[#bucket + 1] = s
			end
		end
	end

	local used: {[any]: boolean} = {}
	local loops: {{any}} = {}
	for _, s0 in ipairs(segs) do
		if not used[s0] then
			local loop, cur = {}, s0
			while cur and not used[cur] do
				used[cur] = true
				loop[#loop + 1] = cur
				local cands = byStart[cur.b[1] .. "," .. cur.b[2]]
				local nxt = nil
				if cands then
					for _, s in ipairs(cands) do
						if not used[s] then nxt = s; break end
					end
				end
				cur = nxt
			end
			if #loop >= 4 then loops[#loops + 1] = loop end
		end
	end
	return loops
end

--------------------------------------------------------------------------
-- STEPS 4-8 — one loop of boundary cells to a clean polygon
--------------------------------------------------------------------------

local function segmentLoop(pts: {P2}, c: any, cls: {string})
	local n = #pts
	local segs: {any} = {}
	local cur = { 1 }
	for i = 2, n do
		-- A run is wall or it is dropoff, never a blend. They are offset by
		-- different amounts, so a segment spanning both has no single correct
		-- push and would have to compromise -- meaning either a wall that is not
		-- stood off from or a ledge that is eroded. Break the run on the class
		-- change and let the two halves be fitted separately.
		if cls[i] ~= cls[cur[1]] then
			local ce, de = fitLine(pts, cur)
			segs[#segs + 1] = { idx = cur, cen = ce, dir = de, class = cls[cur[1]] }
			cur = { i }
			continue
		end
		local trial = table.clone(cur)
		trial[#trial + 1] = i
		local cen, dir = fitLine(pts, trial)
		-- MAXIMUM residual, never the average: an average lets a shallow corner
		-- hide inside a long run, which is precisely the corner worth keeping.
		if maxResidual(pts, trial, cen, dir) > c.fitTol then
			-- The fit just failed. Close the run and start a new one AT this
			-- point, so the two segments share it. A corner is this failure --
			-- a byproduct, not something detected separately against geometry.
			local ce, de = fitLine(pts, cur)
			segs[#segs + 1] = { idx = cur, cen = ce, dir = de, class = cls[cur[1]] }
			cur = { i - 1, i }
		else
			cur = trial
		end
	end
	local ce, de = fitLine(pts, cur)
	segs[#segs + 1] = { idx = cur, cen = ce, dir = de, class = cls[cur[1]] }
	return segs
end

local function spanLength(pts: {P2}, idx: {number}): number
	return len(sub(pts[idx[#idx]], pts[idx[1]]))
end

-- STEP 8. Greedy segmentation is order dependent and the start point on a
-- closed loop is arbitrary, so there is always a spurious corner wherever the
-- walk began: test the first and last segments for collinearity and merge. Then
-- absorb segments too short to be real. Both must happen BEFORE corners are
-- intersected.
local function mergeSegments(pts: {P2}, segs: {any}, c: any)
	local cosLim = math.cos(math.rad(c.collinearDeg))

	-- merge a into b, keeping every point of both, and refuse if the combined
	-- run is no longer straight
	local function refit(a: any, b: any): any?
		-- never across a class change: the merged run would need one push for
		-- what is wall and another for what is ledge
		if a.class ~= b.class then return nil end
		local idx = table.clone(a.idx)
		local seen: {[number]: boolean} = {}
		for _, i in ipairs(idx) do seen[i] = true end
		for _, i in ipairs(b.idx) do
			if not seen[i] then idx[#idx + 1] = i; seen[i] = true end
		end
		local cen, dir = fitLine(pts, idx)
		if maxResidual(pts, idx, cen, dir) > c.fitTol then return nil end
		return { idx = idx, cen = cen, dir = dir, class = a.class }
	end

	-- SEAM. The walk's start point on a closed loop is arbitrary, so a straight
	-- run that happens to straddle it is always split in two -- a spurious
	-- corner at exactly the place nothing happened.
	if #segs >= 2 then
		local first, last = segs[1], segs[#segs]
		if dot(first.dir, last.dir) >= cosLim then
			local m = refit(last, first)
			if m then
				segs[1] = m
				segs[#segs] = nil
			end
		end
	end

	-- Near-collinear and too-short runs, absorbed until nothing changes. Curved
	-- geometry would otherwise shatter into dozens of tiny pieces, and short
	-- pieces are also what makes the corner intersection unstable -- which is
	-- why this runs BEFORE any corner is computed.
	local changed = true
	while changed and #segs > 3 do
		changed = false
		for k = 1, #segs do
			local a, b = segs[k], segs[(k % #segs) + 1]
			if a == b then break end
			-- SIGNED, not abs. Directions are oriented along the walk, so an abs
			-- test calls a 180-degree reversal "collinear" -- and the two long
			-- sides of a thin wall's ring are exactly that, sitting within fitTol
			-- of each other because the wall is thinner than the tolerance. They
			-- merge, the ring falls under three segments, and it collapses to the
			-- raw fallback. That was 165 of 169 holes on SmallMap: every interior
			-- wall came out jagged.
			local nearly = dot(a.dir, b.dir) >= cosLim
			local tiny = spanLength(pts, a.idx) < c.minSegLen
				or spanLength(pts, b.idx) < c.minSegLen
			if nearly or tiny then
				local m = refit(a, b)
				if m then
					segs[k] = m
					table.remove(segs, (k % #segs) + 1)
					changed = true
					break
				end
			end
		end
	end
	return segs
end

--------------------------------------------------------------------------

function Boundary.fromFloor(floorData: any, cfg: Config?)
	local c = merged(cfg)
	local t0 = os.clock()
	local nodes, comps, allY = Boundary.layers(floorData, cfg)

	-- WALL or DROPOFF. Erosion exists to stop an agent clipping a wall or
	-- snagging a corner, and a cliff edge is neither -- an agent walking the lip
	-- of a ledge is fine, and pulling the boundary back from every ledge is what
	-- ate 12.4% of SmallMap's cells. So only wall boundaries are offset.
	--
	-- The test is entirely in the surfel field, no Part consulted: a wall is a
	-- neighbouring column whose surface stands above us, a dropoff is one whose
	-- surface lies below. Nothing at all in that column is open air, which is a
	-- dropoff too. The ambiguous "something at our own height that is not in our
	-- layer" falls to wall, because over-eroding is safe and under-eroding is not.
	local function classify(nk: string, y: number): string
		local ys = allY[nk]
		if not ys then return "drop" end
		local below = false
		for _, h in ipairs(ys) do
			if h > y + c.stepTol then return "wall" end
			if h < y - c.stepTol then below = true end
		end
		return below and "drop" or "wall"
	end

	local stats = {
		cells = #nodes, layers = #comps, regions = 0, holes = 0,
		rawSegments = 0, segments = 0, bevels = 0, unstableCorners = 0,
		rawRings = 0, tinyRegions = 0, wallSegments = 0, dropSegments = 0,
		offsetSum = 0, offsetMin = math.huge, offsetMax = 0,
		severed = 0, severedLayers = 0, annihilated = 0, droppedCells = 0,
		worstResidual = 0,
		heightClamped = 0, heightUnreachable = 0,
	}
	local regions: {any} = {}

	for ci, members in ipairs(comps) do
		-- the layer's own cell mask. A cell can appear twice only where a ramp
		-- spirals over itself inside ONE layer; the mask merges them, which is a
		-- known limit of working in a 2D raster and is recorded, not hidden.
		local member: {[string]: boolean} = {}
		local cells: {any} = {}
		local cellY: {[string]: number} = {}
		local minX, maxX, minZ, maxZ = math.huge, -math.huge, math.huge, -math.huge
		for _, ni in ipairs(members) do
			local nd = nodes[ni]
			local k = nd.ix .. ":" .. nd.iz
			if not member[k] then
				member[k] = true
				cells[#cells + 1] = nd
				cellY[k] = nd.y
			end
			if nd.ix < minX then minX = nd.ix end
			if nd.ix > maxX then maxX = nd.ix end
			if nd.iz < minZ then minZ = nd.iz end
			if nd.iz > maxZ then maxZ = nd.iz end
		end
		if #cells < 3 then continue end

		------------------------------------------------------------------
		-- STEP 2 — distance transform over this layer
		------------------------------------------------------------------
		local ox, oz = minX - 1, minZ - 1
		local w, h = (maxX - minX) + 3, (maxZ - minZ) + 3
		local g = edt2d(w, h, function(i, j)
			return member[(ox + i - 1) .. ":" .. (oz + j - 1)] == true
		end)
		local function distAt(ix: number, iz: number): number
			local i, j = ix - ox + 1, iz - oz + 1
			if i < 1 or i > w or j < 1 or j > h then return 0 end
			return math.sqrt(g[(j - 1) * w + i])
		end

		------------------------------------------------------------------
		-- STEP 3 — contour
		------------------------------------------------------------------
		local loops = traceLoops(member, cells, cellY, c.stepTol)
		local rings: {any} = {}

		-- A ring too small to survive segmentation is still a real obstacle, and
		-- DISCARDING it is the one operation in this module that can invent
		-- walkable ground -- an agent would path straight through a 1x2 pillar.
		-- That breaks the safety property the whole design rests on, so a ring
		-- that cannot be fitted falls back to its raw lattice outline, expanded
		-- away from its own centre by the agent radius. Jagged, but an obstacle
		-- that comes out too big is safe and one that is absent is not.
		local function emitRaw(loop: {any})
			local vs: {P2} = {}
			for _, s in ipairs(loop) do
				local last = vs[#vs]
				if not last or last.x ~= s.a[1] or last.z ~= s.a[2] then
					vs[#vs + 1] = { x = s.a[1], z = s.a[2] }
				end
			end
			if #vs > 1 and vs[1].x == vs[#vs].x and vs[1].z == vs[#vs].z then vs[#vs] = nil end
			if #vs < 3 then return end
			local cx, cz = 0, 0
			for _, p in ipairs(vs) do cx += p.x; cz += p.z end
			cx, cz = cx / #vs, cz / #vs
			for _, p in ipairs(vs) do
				local dx, dz = p.x - cx, p.z - cz
				local m = math.sqrt(dx * dx + dz * dz)
				if m > 1e-6 then
					p.x += dx / m * c.agentRadius
					p.z += dz / m * c.agentRadius
				end
			end
			local a2 = 0
			for i = 1, #vs do
				local p, q = vs[i], vs[(i % #vs) + 1]
				a2 += p.x * q.z - q.x * p.z
			end
			stats.rawRings += 1
			rings[#rings + 1] = { pts = vs, area = a2 * 0.5, raw = true }
		end

		for _, loop in ipairs(loops) do
			-- ordered boundary cell centres. A corner cell contributes two edges;
			-- the duplicate carries no information, so collapse it.
			local pts: {P2} = {}
			local owners: {any} = {}
			local cls: {string} = {}
			local lastCell, lastCls = nil, nil
			for _, s in ipairs(loop) do
				local k = classify(s.nk, s.cell.y)
				-- collapse the duplicate a corner cell contributes, but ONLY while
				-- the class holds: a cell with a wall on one side and a ledge on
				-- the other is two different boundaries and has to stay two points
				if s.cell ~= lastCell or k ~= lastCls then
					pts[#pts + 1] = { x = s.cell.ix + 0.5, z = s.cell.iz + 0.5 }
					owners[#owners + 1] = s.cell
					cls[#cls + 1] = k
					lastCell, lastCls = s.cell, k
				end
			end
			if #pts < 3 then emitRaw(loop); continue end

			-- DESPECKLE THE CLASSIFICATION. The wall/dropoff test is per cell and
			-- reads a neighbouring column, so it picks up grit: a gap between two
			-- parts, a wall whose top surfel is missing, a stair nosing. Left
			-- alone each speck forces the fit to break a run it should have kept,
			-- and since merging refuses to cross a class change the fragments
			-- never recombine -- 854 flips along one layer's 2956 boundary edges,
			-- 38% of them a single edge disagreeing with both its neighbours, and
			-- the segment count nearly doubled. A lone dissenter surrounded by
			-- agreement is grit, so let its neighbours outvote it. Real wall/ledge
			-- transitions are many cells long and unaffected.
			local nC = #cls
			if nC >= 3 then
				for _ = 1, 2 do
					local prev = table.clone(cls)
					for i = 1, nC do
						local a = prev[((i - 2) % nC) + 1]
						local b = prev[(i % nC) + 1]
						if a == b and prev[i] ~= a then cls[i] = a end
					end
				end
			end

			--------------------------------------------------------------
			-- STEP 4 — greedy line fit. THE STAIRCASE DIES HERE.
			--------------------------------------------------------------
			local segs = segmentLoop(pts, c, cls)
			stats.rawSegments += #segs
			segs = mergeSegments(pts, segs, c)
			stats.segments += #segs
			if #segs < 3 then emitRaw(loop); continue end

			--------------------------------------------------------------
			-- STEPS 5 + 6 — bias inward, then offset inward, graded
			--------------------------------------------------------------
			local lines: {any} = {}
			for _, s in ipairs(segs) do
				local out = outwardOf(s.dir)
				-- STEP 5. Translate the line inward until no accepted cell centre
				-- lies outward of it. Without this a fit can sit OUTWARD of the
				-- true wall and eat into the clearance margin, making the safety
				-- guarantee probabilistic instead of exact.
				-- maxD is the LOCAL THICKNESS of the ground behind this line, and it
				-- has to be measured by marching INWARD. Reading D at the boundary
				-- cells themselves is worthless: a boundary cell is 4-adjacent to a
				-- non-walkable cell by construction, so its D is always exactly 1,
				-- which pins the push at clamp(1 - margin, 0, r) forever. The offset
				-- then ignores agentRadius entirely and the grading in step 6 is a
				-- constant. March along the inward normal instead, only as far as
				-- the clamp can still respond to (r + margin); past that the answer
				-- cannot change.
				local reach = math.ceil(c.agentRadius + c.margin)
				local cval, maxD = -math.huge, 0
				for _, i in ipairs(s.idx) do
					local p = pts[i]
					local v = p.x * out.x + p.z * out.z
					if v > cval then cval = v end
					for t = 0, reach do
						local d = distAt(math.floor(p.x - out.x * t), math.floor(p.z - out.z * t))
						if d > maxD then maxD = d end
					end
				end
				local r = maxResidual(pts, s.idx, s.cen, s.dir)
				if r > stats.worstResidual then stats.worstResidual = r end

				-- STEP 6, graded: no threshold, so adjacent polygons on either
				-- side of a narrowing still agree on where the boundary is.
				--
				-- A DROPOFF IS NOT OFFSET AT ALL. The offset exists so an agent
				-- does not clip a wall or snag on a corner; a ledge presents
				-- neither, and walking its lip is legitimate. Standing off from
				-- every ledge as though it were a wall is what removed 12.4% of
				-- SmallMap, and it removes it from exactly the places -- balcony
				-- rims, platform edges, the tops of stairs -- where the ground is
				-- narrowest and most worth keeping.
				local isWall = (s.class or "wall") == "wall"
				local push = isWall and math.clamp(maxD - c.margin, 0, c.agentRadius) or 0
				if isWall then stats.wallSegments += 1 else stats.dropSegments += 1 end
				stats.offsetSum += push
				if push < stats.offsetMin then stats.offsetMin = push end
				if push > stats.offsetMax then stats.offsetMax = push end

				lines[#lines + 1] = {
					n = out, c = cval - push, dir = s.dir,
					anchor = pts[s.idx[#s.idx]], class = s.class or "wall",
				}
			end

			--------------------------------------------------------------
			-- STEP 7 — corners are the intersections of adjacent offset lines
			--
			-- Sub-cell accurate, and sharper than the raster could ever be. This
			-- is also why lines are offset rather than cells eroded up front:
			-- eroding first bevels every convex corner into a small arc, which
			-- the segmenter then reads as two or three short segments instead of
			-- one corner -- more polygons, the opposite of the goal.
			--------------------------------------------------------------
			local verts: {P2} = {}
			local nL = #lines
			for i = 1, nL do
				local l1, l2 = lines[i], lines[(i % nL) + 1]
				local det = l1.n.x * l2.n.z - l1.n.z * l2.n.x
				local anchor = l1.anchor
				if math.abs(det) < 1e-6 then
					-- near-parallel: fall back to the foot of the anchor on l1
					stats.unstableCorners += 1
					local s = l1.c - (anchor.x * l1.n.x + anchor.z * l1.n.z)
					verts[#verts + 1] = { x = anchor.x + l1.n.x * s, z = anchor.z + l1.n.z * s }
				else
					local px = (l1.c * l2.n.z - l2.c * l1.n.z) / det
					local pz = (l1.n.x * l2.c - l2.n.x * l1.c) / det
					local d = len(sub({ x = px, z = pz }, anchor))
					if d > c.miterLimit then
						-- MITER LIMIT. An acute corner throws the intersection
						-- arbitrarily far out; bevel across it instead.
						stats.bevels += 1
						local s1 = l1.c - (anchor.x * l1.n.x + anchor.z * l1.n.z)
						verts[#verts + 1] = { x = anchor.x + l1.n.x * s1, z = anchor.z + l1.n.z * s1 }
						local s2 = l2.c - (anchor.x * l2.n.x + anchor.z * l2.n.z)
						verts[#verts + 1] = { x = anchor.x + l2.n.x * s2, z = anchor.z + l2.n.z * s2 }
					else
						verts[#verts + 1] = { x = px, z = pz }
					end
				end
			end
			if #verts < 3 then emitRaw(loop); continue end

			-- signed area in XZ decides outer vs hole. Magnitude picks the outer
			-- ring, never the sign: this project has been bitten by a handedness
			-- assumption before.
			local a2 = 0
			for i = 1, #verts do
				local p, q = verts[i], verts[(i % #verts) + 1]
				a2 += p.x * q.z - q.x * p.z
			end
			rings[#rings + 1] = { pts = verts, area = a2 * 0.5 }
		end

		if #rings == 0 then continue end
		local outer = rings[1]
		for _, r in ipairs(rings) do
			if math.abs(r.area) > math.abs(outer.area) then outer = r end
		end
		-- The raw fallback expands a ring away from its own centre, which is the
		-- safe direction for a HOLE (an obstacle that comes out too big costs
		-- nothing) and exactly the wrong one for an OUTER ring, where it would
		-- hand back more walkable ground than exists. A ring that could not be
		-- segmented is at most a couple of studs across, so as an outer ring it
		-- describes a ledge no agent of radius r could stand on. Drop the region
		-- and count it: removing ground is erosion, and erosion is always safe.
		if outer.raw then
			stats.tinyRegions += 1
			continue
		end
		local holes = {}
		for _, r in ipairs(rings) do if r ~= outer then holes[#holes + 1] = r end end

		-- height for a vertex: the nearest cell of THIS layer. Vertical accuracy
		-- is already exact per cell, so nothing is interpolated across a step.
		-- A vertex is offset off its own cells by up to the agent radius, and a raw
		-- ring is expanded further still, so it can easily land on a cell this
		-- layer does not own. The old last resort was cells[1].y -- an arbitrary
		-- member of the layer, which on a layer spanning 15 studs of ramp is an
		-- arbitrary height. That is what drew the vertical fences at the wall
		-- bases: one vertex of a flat ring teleporting to whatever height cell 1
		-- happened to sit at. Search wider, then fall back to the genuinely
		-- nearest cell rather than the first one -- and, per `toWorld` below,
		-- never further vertically than a walkable slope could have carried the
		-- ring there.
		-- Nearest cell of this layer whose height falls inside [lo, hi]. Returns
		-- nil when the layer has no such cell anywhere, so the caller can decide
		-- what to do rather than being handed a wrong number.
		local function nearestCell(p: P2, lo: number, hi: number): (number?, number)
			local bx, bz = math.floor(p.x), math.floor(p.z)
			for rad = 0, 8 do
				local best, bestD = nil, math.huge
				for dx = -rad, rad do
					for dz = -rad, rad do
						if math.max(math.abs(dx), math.abs(dz)) == rad then
							local k = (bx + dx) .. ":" .. (bz + dz)
							local y = cellY[k]
							if y and y >= lo and y <= hi then
								local d = (bx + dx + 0.5 - p.x) ^ 2 + (bz + dz + 0.5 - p.z) ^ 2
								if d < bestD then best, bestD = y, d end
							end
						end
					end
				end
				if best then return best, bestD end
			end
			local best, bestD = nil, math.huge
			for _, cell in ipairs(cells) do
				if cell.y >= lo and cell.y <= hi then
					local d = (cell.ix + 0.5 - p.x) ^ 2 + (cell.iz + 0.5 - p.z) ^ 2
					if d < bestD then best, bestD = cell.y, d end
				end
			end
			return best, bestD
		end

		-- VERTICAL CONTINUITY, and it is what killed the picket fences.
		--
		-- Nearest-cell-in-XZ is the right height only while the layer is single
		-- valued near the vertex. A layer legitimately spans a whole building --
		-- ground, stairs and roof are one connected, x:z-injective component --
		-- so a vertex sitting on the roof's edge has ground cells 1 stud away in
		-- XZ and 10 studs below. Nearest-in-XZ picked those, and the ring
		-- teleported down the wall and back: on SmallMap 58 of 102 holes had two
		-- CONSECUTIVE vertices, 2-4 studs apart on the ground, 10 studs apart in
		-- height. Drawn, that is a fence of vertical bars down every facade.
		--
		-- A ring is a loop of walkable boundary, so the constraint is just the
		-- walkable slope: over dxz studs of ground the height may move by at most
		-- stepTol per stud, which is the same 65-degree limit layering uses. Walk
		-- the loop from the vertex whose unconstrained answer is most trustworthy
		-- -- the one sitting closest to a cell of its own layer -- and let each
		-- vertex pick the nearest cell that its predecessor's height can reach.
		--
		-- The constraint is never allowed to FAIL a vertex: if no cell in the
		-- layer is reachable, the unconstrained nearest is still better than
		-- nothing, and the fallback is counted rather than hidden.
		local function toWorld(ring: any): {Vector3}
			local pts = ring.pts
			local n = #pts
			local free, freeD = {}, {}
			for i = 1, n do
				local y, d = nearestCell(pts[i], -math.huge, math.huge)
				free[i], freeD[i] = y or cells[1].y, d
			end
			local seed, seedD = 1, math.huge
			for i = 1, n do
				if freeD[i] < seedD then seed, seedD = i, freeD[i] end
			end

			local out = table.create(n)
			out[seed] = Vector3.new(pts[seed].x, free[seed], pts[seed].z)
			local prevY = free[seed]
			for step = 1, n - 1 do
				local i = ((seed - 1 + step) % n) + 1
				local prev = ((seed - 2 + step) % n) + 1
				local rise = c.stepTol * math.max(1, len(sub(pts[i], pts[prev])))
				local y = nearestCell(pts[i], prevY - rise, prevY + rise)
				if y then
					if y ~= free[i] then stats.heightClamped += 1 end
				else
					y = free[i]
					stats.heightUnreachable += 1
				end
				out[i] = Vector3.new(pts[i].x, y, pts[i].z)
				prevY = y
			end
			return out
		end

		local region = {
			layer = ci,
			verts = toWorld(outer),
			holes = {},
			area = math.abs(outer.area),
			cells = #cells,
		}
		for _, hr in ipairs(holes) do
			region.holes[#region.holes + 1] = { verts = toWorld(hr), area = math.abs(hr.area) }
		end
		regions[#regions + 1] = region
		stats.regions += 1
		stats.holes += #region.holes

		------------------------------------------------------------------
		-- STEP 9 — severance check
		--
		-- Mandatory, not optional: this is the ONLY detector for the offset
		-- having pinched a corridor shut and cut the map in half. Union-find
		-- over surfel adjacency before and after; anything that loses
		-- connectivity is reported. Reported, never silently repaired -- a
		-- severed layer is a tuning failure and must be visible as one.
		------------------------------------------------------------------
		local function inRing(ring: {P2}, x: number, z: number): boolean
			local inside = false
			local n = #ring
			local j = n
			for i = 1, n do
				local a, b = ring[i], ring[j]
				if (a.z > z) ~= (b.z > z) then
					local t = (z - a.z) / (b.z - a.z)
					if x < a.x + t * (b.x - a.x) then inside = not inside end
				end
				j = i
			end
			return inside
		end
		local kept: {[string]: boolean} = {}
		local keptList: {any} = {}
		for _, cell in ipairs(cells) do
			local x, z = cell.ix + 0.5, cell.iz + 0.5
			local ok = inRing(outer.pts, x, z)
			if ok then
				for _, hr in ipairs(holes) do
					if inRing(hr.pts, x, z) then ok = false; break end
				end
			end
			if ok then
				kept[cell.ix .. ":" .. cell.iz] = true
				keptList[#keptList + 1] = cell
			end
		end
		stats.droppedCells += (#cells - #keptList)
		local seen: {[string]: boolean} = {}
		local pieces = 0
		local D4 = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } }
		for _, cell in ipairs(keptList) do
			local k0 = cell.ix .. ":" .. cell.iz
			if not seen[k0] then
				seen[k0] = true
				pieces += 1
				local queue, head = { cell }, 1
				while head <= #queue do
					local cur = queue[head]; head += 1
					for _, d in ipairs(D4) do
						local nk = (cur.ix + d[1]) .. ":" .. (cur.iz + d[2])
						if kept[nk] and not seen[nk] then
							seen[nk] = true
							queue[#queue + 1] = { ix = cur.ix + d[1], iz = cur.iz + d[2] }
						end
					end
				end
			end
		end
		-- Counting the pieces of what survives cannot see a layer that did not
		-- survive at all: annihilation is zero pieces, and `pieces > 1` is false
		-- for zero exactly as it is for one. A corridor the offset ate whole is
		-- the most severe version of the failure this check exists for, so it has
		-- to be the loudest case, not the one that slips past. A 4-stud corridor
		-- does this on the default tuning: maxD is 2 at the centreline, so both
		-- sides push the full 1.5 and the two offset lines land on each other.
		if #keptList == 0 then
			stats.severedLayers += 1
			stats.annihilated += 1
			region.severedInto = 0
		elseif pieces > 1 then
			stats.severedLayers += 1
			stats.severed += (pieces - 1)
			region.severedInto = pieces
		end
	end

	----------------------------------------------------------------------
	-- STEP 10 — inter-layer links
	--
	-- Layers are 2.5D rasters and a raster cannot express "upstairs". Splitting
	-- by connectivity AND x:z injectivity is what makes each raster meaningful,
	-- but it also cuts the ramp free of the floor it climbs from, so the output
	-- is a pile of disconnected islands: correct surfaces, no way between them.
	--
	-- The cut is recoverable exactly, because it is the same predicate that made
	-- it. Two cells 4-adjacent in XZ and within a step of each other are
	-- walkable neighbours; if they landed in different layers, that is precisely
	-- where a layer boundary was drawn through traversable ground -- a ramp foot,
	-- a stair head, a balcony meeting its walkway. A drop bigger than stepTol is
	-- NOT a link, which is why cliffs stay cliffs.
	--
	-- Crossings are clustered into runs so a pair of layers touching in two
	-- places yields two links rather than one averaged nonsense in between.
	----------------------------------------------------------------------
	local atAll: {[string]: {number}} = {}
	for i, nd in ipairs(nodes) do
		local k = nd.ix .. ":" .. nd.iz
		local b = atAll[k]
		if not b then b = {}; atAll[k] = b end
		b[#b + 1] = i
	end

	local D4L = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } }
	local crossings: {[string]: {any}} = {}
	for i, nd in ipairs(nodes) do
		for _, d in ipairs(D4L) do
			local b = atAll[(nd.ix + d[1]) .. ":" .. (nd.iz + d[2])]
			if b then
				for _, ni in ipairs(b) do
					local nb = nodes[ni]
					-- one direction only, or every crossing is found twice
					if nb.comp > nd.comp and math.abs(nb.y - nd.y) <= c.stepTol then
						local pk = nd.comp .. ":" .. nb.comp
						local list = crossings[pk]
						if not list then list = {}; crossings[pk] = list end
						list[#list + 1] = {
							ix = nd.ix, iz = nd.iz,
							from = Vector3.new(nd.ix + 0.5, nd.y, nd.iz + 0.5),
							to = Vector3.new(nb.ix + 0.5, nb.y, nb.iz + 0.5),
						}
					end
				end
			end
		end
	end

	local links: {any} = {}
	for pk, list in pairs(crossings) do
		local a, bId = pk:match("(%d+):(%d+)")
		-- cluster by adjacency of the crossing cells, so two separate staircases
		-- between the same pair of layers stay two links
		local byCell: {[string]: {number}} = {}
		for idx, cr in ipairs(list) do
			local k = cr.ix .. ":" .. cr.iz
			local bb = byCell[k]
			if not bb then bb = {}; byCell[k] = bb end
			bb[#bb + 1] = idx
		end
		local taken: {[number]: boolean} = {}
		for idx = 1, #list do
			if not taken[idx] then
				local group = { idx }
				taken[idx] = true
				local queue, head = { idx }, 1
				while head <= #queue do
					local cur = list[queue[head]]; head += 1
					for _, d in ipairs(D4L) do
						local bb = byCell[(cur.ix + d[1]) .. ":" .. (cur.iz + d[2])]
						if bb then
							for _, j in ipairs(bb) do
								if not taken[j] then
									taken[j] = true
									group[#group + 1] = j
									queue[#queue + 1] = j
								end
							end
						end
					end
				end
				local fx, fy, fz, tx, ty, tz = 0, 0, 0, 0, 0, 0
				for _, j in ipairs(group) do
					local cr = list[j]
					fx += cr.from.X; fy += cr.from.Y; fz += cr.from.Z
					tx += cr.to.X;   ty += cr.to.Y;   tz += cr.to.Z
				end
				local n = #group
				links[#links + 1] = {
					a = tonumber(a), b = tonumber(bId),
					cells = n,
					-- how wide the opening is, in studs of shared frontage
					width = n,
					from = Vector3.new(fx / n, fy / n, fz / n),
					to = Vector3.new(tx / n, ty / n, tz / n),
				}
			end
		end
	end
	stats.links = #links
	local linkedLayers: {[number]: boolean} = {}
	for _, l in ipairs(links) do linkedLayers[l.a] = true; linkedLayers[l.b] = true end
	local nLinked = 0
	for _ in pairs(linkedLayers) do nLinked += 1 end
	stats.linkedLayers = nLinked
	stats.isolatedLayers = #comps - nLinked

	if stats.offsetMin == math.huge then stats.offsetMin = 0 end
	stats.seconds = os.clock() - t0
	return { regions = regions, stats = stats, config = c, links = links }
end

--------------------------------------------------------------------------

function Boundary.visualize(res: any, parent: Instance?)
	local root = parent or workspace
	local old = root:FindFirstChild("NavGen_Boundary")
	if old then old:Destroy() end
	local folder = Instance.new("Folder")
	folder.Name = "NavGen_Boundary"
	folder.Parent = root

	local function bar(a: Vector3, b: Vector3, colour: Color3, thick: number, into: Instance)
		local d = b - a
		if d.Magnitude < 1e-4 then return end
		local p = Instance.new("Part")
		p.Anchored = true; p.CanCollide = false; p.CanQuery = false; p.CanTouch = false
		p.Material = Enum.Material.Neon
		p.Color = colour
		p.Size = Vector3.new(thick, thick, d.Magnitude)
		p.CFrame = CFrame.lookAt(a + d * 0.5, b) + Vector3.new(0, 0.08, 0)
		p.Parent = into
	end

	for ri, r in ipairs(res.regions) do
		local sub = Instance.new("Folder")
		sub.Name = string.format("R%d_L%d_c%d_h%d%s", ri, r.layer, r.cells, #r.holes,
			r.severedInto and ("_SEVERED" .. r.severedInto) or "")
		sub.Parent = folder
		local v = r.verts
		for i = 1, #v do
			bar(v[i], v[(i % #v) + 1], Color3.fromRGB(255, 80, 80), 0.18, sub)
		end
		for _, hr in ipairs(r.holes) do
			local hv = hr.verts
			for i = 1, #hv do
				bar(hv[i], hv[(i % #hv) + 1], Color3.fromRGB(80, 200, 255), 0.18, sub)
			end
		end
	end

	-- Links, in green: where one layer hands over to another. If a staircase
	-- has no green bar at its foot or its head, the bake has produced two
	-- surfaces with no way between them, and that is visible at a glance.
	if res.links and #res.links > 0 then
		local lf = Instance.new("Folder")
		lf.Name = "Links"
		lf.Parent = folder
		for _, l in ipairs(res.links) do
			local seg = Instance.new("Folder")
			seg.Name = string.format("L%d_%d_w%d", l.a, l.b, l.width)
			seg.Parent = lf
			-- from -> to is one cell step, so drawing that is a bar with no
			-- length. What matters is how WIDE the opening is: draw it across
			-- the frontage, perpendicular to the direction of travel.
			local mid = (l.from + l.to) * 0.5 + Vector3.new(0, 0.5, 0)
			local step = l.to - l.from
			local flat = Vector3.new(step.X, 0, step.Z)
			local across = (flat.Magnitude > 1e-3)
				and Vector3.new(-flat.Unit.Z, 0, flat.Unit.X)
				or Vector3.new(1, 0, 0)
			local half = across * (l.width * 0.5)
			bar(mid - half, mid + half, Color3.fromRGB(90, 255, 120), 0.5, seg)
		end
	end
	return folder
end

return Boundary
